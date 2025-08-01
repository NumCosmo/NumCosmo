#!/usr/bin/env python3
"""Generate stubs for numcosmo_py."""

# pylint: disable=too-many-lines
# Based on https://github.com/PyCQA/astroid/blob/main/astroid/brain/brain_gi.py
# Licensed under the LGPL: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html
# For details: https://github.com/PyCQA/astroid/blob/main/LICENSE
# Updated for numcosmo_py by Sandro Dias Pinto Vitenti <vitenti@uel.br> in July 2025

from __future__ import annotations

from typing import Any, Callable, Optional, Tuple, Type, Union, cast

import argparse
import importlib
import inspect
import itertools
import pprint
import re
import textwrap
from types import ModuleType

import gi
import gi._gi as GI
from numcosmo_py import parse_pyi

gi.require_version("GIRepository", "3.0")

# pyright: reportMissingModuleSource=false
# pylint:disable-next=wrong-import-position,unused-import,wrong-import-order
from gi.repository import GIRepository  # noqa: E402

# pyright: reportMissingModuleSource=false
# pylint:disable-next=wrong-import-position,unused-import,wrong-import-order
from gi.repository import GObject  # noqa: E402

_identifier_re = r"^[A-Za-z_]\w*$"

ObjectT = Union[ModuleType, Type[Any]]  # pylint: disable=invalid-name
dunder_list = ["__getitem__", "__setitem__"]


def _object_get_props(
    repo: GIRepository.Repository, obj: GI.ObjectInfo
) -> Tuple[list[GIRepository.BaseInfo], list[GIRepository.BaseInfo]]:
    parents: list[GI.ObjectInfo] = []
    parent: Optional[GI.ObjectInfo] = obj.get_parent()
    while parent:
        parents.append(parent)
        parent = parent.get_parent()

    interfaces: list[GI.InterfaceInfo] = list(obj.get_interfaces())

    subclasses: list[GI.ObjectInfo | GI.InterfaceInfo] = parents + interfaces

    props: list[GI.PropertyInfo] = list(obj.get_properties())
    for s in subclasses:
        props.extend(s.get_properties())

    readable_props: list[GIRepository.BaseInfo] = []
    writable_props: list[GIRepository.BaseInfo] = []
    for prop in props:
        namespace = prop.get_namespace()
        container = prop.get_container()
        class_info = repo.find_by_name(namespace, container.get_name())
        if class_info is None:
            raise RuntimeError(f"Unable to find {namespace}.{container}")

        assert isinstance(class_info, GIRepository.ObjectInfo)
        # pylint: disable-next=no-member
        if class_info.get_g_type().is_a(GObject.Object.__gtype__):  # type: ignore
            n_props = class_info.get_n_properties()
            for i in range(n_props):
                p = class_info.get_property(i)
                if p.get_name() == prop.get_name():
                    flags = p.get_flags()
                    if flags & GObject.ParamFlags.READABLE:
                        readable_props.append(p)
                    if flags & GObject.ParamFlags.WRITABLE:
                        writable_props.append(p)

        if class_info.get_g_type().is_interface():
            n_props = class_info.get_n_properties()
            for i in range(n_props):
                p = class_info.get_property(i)
                if p.get_name() == prop.get_name():
                    flags = p.get_flags()
                    if flags & GObject.ParamFlags.READABLE:
                        readable_props.append(p)
                    if flags & GObject.ParamFlags.WRITABLE:
                        writable_props.append(p)

    return (readable_props, writable_props)


def _callable_get_arguments(
    ltype: GI.CallbackInfo,
    current_namespace: str,
    needed_namespaces: set[str],
    can_default: bool = False,
) -> Tuple[list[str], list[str], list[str]]:
    function_args = ltype.get_arguments()
    accept_optional_args = False
    optional_args_name = ""
    dict_names: dict[int, str] = {}
    dict_args: dict[int, GI.ArgInfo] = {}
    str_args: list[str] = []
    dict_return_args: dict[int, str] = {}
    skip: list[int] = []

    # Filter out array length arguments for return ltype
    ret_type = ltype.get_return_type()
    if ret_type.get_array_length_index() >= 0:
        skip.append(ret_type.get_array_length_index())

    for i, arg in enumerate(function_args):
        if i in skip:
            continue

        if arg.get_closure_index() >= 0:
            accept_optional_args = True
            optional_args_name = function_args[arg.get_closure_index()].get_name()
            skip.append(arg.get_closure_index())
            skip.append(arg.get_destroy_index())

        # Filter out array length args
        arg_type = arg.get_type_info()
        len_arg: int = arg_type.get_array_length_index()
        if len_arg >= 0:
            skip.append(len_arg)
            if len_arg < i:
                dict_names.pop(len_arg, None)
                dict_args.pop(len_arg, None)
                dict_return_args.pop(len_arg, None)

        # Need to check because user_data can be the first arg
        if arg.get_closure_index() != i and arg.get_destroy_index() != i:
            direction = arg.get_direction()
            if direction in (GI.Direction.OUT, GI.Direction.INOUT):
                t = _type_to_python(
                    arg.get_type_info(), current_namespace, needed_namespaces, True
                )

                dict_return_args[i] = t
            elif direction in (GI.Direction.IN, GI.Direction.INOUT):
                dict_names[i] = arg.get_name()
                dict_args[i] = arg

    # Traverse args in reverse to check for optional args
    local_args = list(dict_args.values())
    for a in reversed(local_args):
        t = _type_to_python(
            a.get_type_info(),
            current_namespace,
            needed_namespaces,
            False,
            a.get_closure_index() >= 0,  # True if function admits variable arguments
        )

        if a.may_be_null() and t != "None":
            if can_default:
                str_args.append(f"typing.Optional[{t}] = None")
            else:
                str_args.append(f"typing.Optional[{t}]")
        else:
            can_default = False
            str_args.append(t)

    str_args = list(reversed(str_args))
    names = list(dict_names.values())

    # We need this info to filter out None as return arg for methods
    # that process Gio.AsyncResult. In python this method raises always.
    is_async_res = "Gio.AsyncResult" in str_args

    if accept_optional_args:
        names.append(f"*{optional_args_name}")
        str_args.append("typing.Any")

    return_type = _type_to_python(
        ltype.get_return_type(), current_namespace, needed_namespaces, True
    )
    if ltype.may_return_null() and return_type != "None" and not is_async_res:
        return_type = f"typing.Optional[{return_type}]"

    return_args = list(dict_return_args.values())
    if return_type != "None" or len(return_args) == 0:
        return_args.insert(0, return_type)

    return (names, str_args, return_args)


class TypeInfo:
    """This struct tries to emulate gi.TypeInfo."""

    def __init__(
        self,
        obj: Any,
        get_tag: Callable[[TypeInfo], int],
        get_param_type: Callable[[TypeInfo, int], TypeInfo],
        get_interface: Callable[[TypeInfo], TypeInfo],
    ):
        """Construct."""
        self.obj = obj
        self._get_tag = get_tag
        self._get_param_type = get_param_type
        self._get_interface = get_interface

    def get_tag(self) -> int:
        """Get tag."""
        return self._get_tag(self.obj)

    def get_param_type(self, n: int) -> TypeInfo:
        """Get param type."""
        ltype = self._get_param_type(self.obj, n)
        return TypeInfo(ltype, self._get_tag, self._get_param_type, self._get_interface)

    def get_interface(self) -> TypeInfo:
        """Get interface."""
        ltype = self._get_interface(self.obj)
        return TypeInfo(ltype, self._get_tag, self._get_param_type, self._get_interface)

    def get_name(self) -> str:
        """Get name."""
        return self.obj.get_name()

    def get_namespace(self) -> str:
        """Get namespace."""
        return self.obj.get_namespace()

    def get_type(self) -> GIRepository.TypeInfo:
        """Get type."""
        return self.obj.get_type()


def _type_to_python(
    ltype: GIRepository.TypeInfo,
    current_namespace: str,
    needed_namespaces: set[str],
    out_arg: bool = False,
    varargs: bool = False,
) -> str:
    tag = ltype.get_tag()
    tags = GI.TypeTag

    if tag == tags.ARRAY:
        array_type = ltype.get_param_type(0)
        assert array_type is not None
        t = _type_to_python(array_type, current_namespace, needed_namespaces)
        if out_arg:
            # As output argument array of ltype uint8 are returned as bytes
            if array_type.get_tag() == GI.TypeTag.UINT8:
                return "bytes"
            return f"list[{t}]"

        # As input arguments array can be generated by any sequence
        if t == "float":
            return "typing.Sequence[float] | npt.NDArray[np.float64]"
        return f"typing.Sequence[{t}]"

    if tag in (tags.GLIST, tags.GSLIST):
        array_type = ltype.get_param_type(0)
        assert array_type is not None
        t = _type_to_python(array_type, current_namespace, needed_namespaces)
        return f"list[{t}]"

    if tag == tags.BOOLEAN:
        return "bool"

    if tag in (tags.DOUBLE, tags.FLOAT):
        return "float"

    if tag == tags.ERROR:
        if current_namespace == "GLib":
            return "Error"
        needed_namespaces.add("GLib")
        return "GLib.Error"

    if tag == tags.GHASH:
        key_type = ltype.get_param_type(0)
        value_type = ltype.get_param_type(1)
        assert key_type is not None
        assert value_type is not None

        kt = _type_to_python(key_type, current_namespace, needed_namespaces)
        vt = _type_to_python(value_type, current_namespace, needed_namespaces)
        return f"dict[{kt}, {vt}]"

    if tag in (tags.FILENAME, tags.UTF8, tags.UNICHAR):
        return "str"

    if tag == tags.GTYPE:
        return "typing.Type[typing.Any]"

    if tag in (
        tags.INT8,
        tags.INT16,
        tags.INT32,
        tags.INT64,
        tags.UINT8,
        tags.UINT16,
        tags.UINT32,
        tags.UINT64,
    ):
        return "int"

    if tag == tags.INTERFACE:
        interface = ltype.get_interface()
        if isinstance(interface, GI.CallbackInfo):
            (names, args, return_args) = _callable_get_arguments(
                interface, current_namespace, needed_namespaces
            )

            return_type = ""
            if len(return_args) == 1:
                return_type = return_args[0]
            else:
                return_type = f"typing.Tuple[{', '.join(return_args)}]"

            # FIXME, how to express Callable with variable arguments?
            if (len(names) > 0 and names[-1].startswith("*")) or varargs:
                return f"typing.Callable[..., {return_type}]"
            return f"typing.Callable[[{', '.join(args)}], {return_type}]"
        else:
            assert isinstance(interface, (GIRepository.BaseInfo, GI.BaseInfo))
            namespace = interface.get_namespace()
            name = interface.get_name()
            assert name is not None

            if not re.match(_identifier_re, name):
                raise ValueError(f"Invalid interface name: {name}")

            if namespace == "GObject" and name == "Value":
                return "typing.Any"

            if namespace == "GObject" and name == "Closure":
                return "typing.Callable[..., typing.Any]"

            if namespace == "cairo" and name == "Context" and not out_arg:
                return "cairo.Context[_SomeSurface]"

            if current_namespace == namespace:
                return f"{name}"
            else:
                needed_namespaces.add(namespace)
                return f"{namespace}.{name}"

    if tag == tags.VOID:
        return "None"

    raise ValueError("TODO")


def _build(
    repo: GIRepository.Repository,
    parent: ObjectT,
    namespace: str,
    local_overrides: dict[str, str],
) -> str:
    ns: set[str] = set()
    ret = _gi_build_stub(
        repo, parent, namespace, dir(parent), ns, local_overrides, None, ""
    )

    typevars: list[str] = [
        'T = typing.TypeVar("T")',
    ]

    if namespace == "Gtk":
        typevars.append(
            """CellRendererT = typing.TypeVar(
    "CellRendererT",
    CellRendererCombo,
    CellRendererPixbuf,
    CellRendererProgress,
    CellRendererSpin,
    CellRendererSpinner,
    CellRendererText,
    CellRendererToggle,
)"""
        )

    imports: list[str] = []
    if "cairo" in ns:
        imports = ["import cairo"]
        typevars.append(
            '_SomeSurface = typing.TypeVar("_SomeSurface", bound=cairo.Surface)'
        )
        ns.remove("cairo")

    imports += [f"from gi.repository import {n}" for n in sorted(ns)]

    return (
        "import typing"
        + "\n\n"
        + "import numpy.typing as npt"
        + "\n\n"
        + "\n".join(imports)
        + "\n"
        + "\n".join(typevars)
        + "\n\n"
        + ret
    )


def _generate_full_name(prefix: str, name: str) -> str:
    full_name = name
    if len(prefix) > 0:
        full_name = f"{prefix}.{name}"
    return full_name


def _build_function_info(
    current_namespace: str,
    name: str,
    function: GI.FunctionInfo | GI.VFuncInfo,
    in_class: Optional[Any],
    needed_namespaces: set[str],
    return_signature: Optional[str] = None,
    comment: Optional[str] = None,
) -> str:
    constructor: bool = False
    method: bool = isinstance(function, GI.VFuncInfo)
    static: bool = False

    # Flags
    function_flags = function.get_flags()
    if function_flags & GI.FunctionInfoFlags.IS_CONSTRUCTOR:
        constructor = True

    if function_flags & GI.FunctionInfoFlags.IS_METHOD:
        method = True

    if in_class and not method and not constructor:
        static = True

    # Arguments
    (names, args, return_args) = _callable_get_arguments(
        function, current_namespace, needed_namespaces, True
    )
    args_types = [f"{name}: {args[i]}" for (i, name) in enumerate(names)]

    # Return ltype
    if return_signature:
        return_type = return_signature
    elif len(return_args) > 1:
        return_type = f"typing.Tuple[{', '.join(return_args)}]"
    else:
        return_type = f"{return_args[0]}"

    # Generate string
    prepend = ""
    if constructor:
        args_types.insert(0, "cls")
        prepend = "@classmethod\n"
        # Override return value, for example Gtk.Button.new returns a Gtk.Widget instead
        # of Gtk.Button
        rt = function.get_container().get_name()
        if return_type != f"typing.Optional[{rt}]":
            return_type = rt
    elif method:
        args_types.insert(0, "self")
    elif static:
        prepend = "@staticmethod\n"

    if comment:
        return (
            f"{prepend}def {name}({', '.join(str(a) for a in args_types)})"
            f" -> {return_type}: ... # {comment}\n"
        )
    else:
        return (
            f"{prepend}def {name}({', '.join(str(a) for a in args_types)})"
            f" -> {return_type}: ...\n"
        )


def _wrapped_strip_boolean_result(
    current_namespace: str,
    name: str,
    function: Any,
    in_class: Optional[Any],
    needed_namespaces: set[str],
) -> str:
    real_function = function.__wrapped__
    fail_ret = inspect.getclosurevars(function).nonlocals.get("fail_ret")

    (_, _, return_args) = _callable_get_arguments(
        real_function, current_namespace, needed_namespaces
    )
    return_args = return_args[1:]  # Strip first return value

    if len(return_args) > 1:
        return_signature = f"typing.Tuple[{', '.join(return_args)}]"
    else:
        return_signature = f"{return_args[0]}"

    if fail_ret is None:
        return_signature = f"typing.Optional[{return_signature}]"
    else:
        if isinstance(fail_ret, tuple):
            if len(fail_ret) > 0:
                return_signature = (
                    f"({return_signature} | "
                    f"typing.Tuple{str(fail_ret).replace('(', '[').replace(')', ']')})"
                )
            else:
                return_signature = f"({return_signature} | typing.Tuple[()])"
        else:
            return_signature = f"({return_signature} | typing.Literal[{fail_ret}])"

    return _build_function_info(
        current_namespace,
        name,
        real_function,
        in_class,
        needed_namespaces,
        return_signature,
        "CHECK Wrapped function",
    )


def _build_function(
    current_namespace: str,
    name: str,
    function: Any,
    in_class: Optional[Any],
    needed_namespaces: set[str],
) -> str:
    if name.startswith("_") and name not in dunder_list:
        return ""

    if hasattr(function, "__wrapped__"):
        if "strip_boolean_result" in str(function):
            return _wrapped_strip_boolean_result(
                current_namespace, name, function, in_class, needed_namespaces
            )

    if isinstance(function, (GI.FunctionInfo, GI.VFuncInfo)):
        return _build_function_info(
            current_namespace, name, function, in_class, needed_namespaces
        )
    if hasattr(function, "__func__") and isinstance(
        function.__func__, (GI.FunctionInfo, GI.VFuncInfo)
    ):
        return _build_function_info(
            current_namespace, name, function.__func__, in_class, needed_namespaces
        )

    signature = str(inspect.signature(function))
    definition = f"def {name}{signature}: ... # FIXME Function\n"

    return definition


def _check_override(
    prefix: str, name: str, local_overrides: dict[str, str]
) -> Optional[str]:
    full_name = _generate_full_name(prefix, name)
    if full_name in local_overrides:
        return "# override\n" + local_overrides[full_name]
    return None


def _gi_build_stub(
    repo: GIRepository.Repository,
    parent: ObjectT,
    current_namespace: str,
    children: list[str],
    needed_namespaces: set[str],
    local_overrides: dict[str, str],
    in_class: Optional[Any],
    prefix_name: str,
) -> str:
    """Build stubs for a GI module.

    Inspect the passed module recursively and build stubs for functions, classes, etc.
    """
    classes: dict[str, Type[Any]] = {}
    functions: dict[str, Callable[..., Any]] = {}
    constants: dict[str, Any] = {}
    flags: dict[str, Type[Any]] = {}
    enums: dict[str, Type[Any]] = {}

    ret = ""

    for name in children:
        if name not in dunder_list:
            if name.startswith("__"):
                continue
            # Check if this is a valid name in python
            if not re.match(_identifier_re, name):
                continue

        try:
            obj = getattr(parent, name)
        except AttributeError:
            continue

        if inspect.isclass(obj):
            if GObject.GFlags in obj.__mro__:
                flags[name] = obj
            elif GObject.GEnum in obj.__mro__:
                enums[name] = obj
            else:
                classes[name] = obj
        elif inspect.isfunction(obj) or inspect.isbuiltin(obj):
            functions[name] = obj
        elif inspect.ismethod(obj) or inspect.ismethoddescriptor(obj):
            functions[name] = obj
        elif callable(obj):
            # Fall back to a function for anything callable
            functions[name] = obj
        elif in_class:
            # Check if obj was already processed
            if hasattr(in_class, "__info__"):
                obj_info = in_class.__info__
                if isinstance(obj_info, (GI.StructInfo, GI.ObjectInfo)):
                    if not (name in [f.get_name() for f in obj_info.get_fields()]):
                        constants[name] = obj
                else:
                    constants[name] = obj
            else:
                constants[name] = obj
        else:
            # Assume everything else is some manner of constant
            constants[name] = obj

    # Constants
    for name in sorted(constants):
        if name[0].isdigit():
            # GDK has some busted constant names like
            # Gdk.EventType.2BUTTON_PRESS
            continue

        override = _check_override(prefix_name, name, local_overrides)
        if override:
            ret += override + "\n"
            continue

        val = constants[name]

        if str(val).startswith(("<flags", "<enum")):
            val = val.real

        if isinstance(val, str):
            ret += f'{name}: {val.__class__.__name__} = r"{val}"\n'
        elif isinstance(val, (bool, float, int)):
            ret += f"{name}: {val.__class__.__name__} = {val}\n"
        elif val.__class__.__name__ == "Atom":
            ret += f"{name}: {val.__class__.__name__} = ...\n"
        else:
            ret += f"{name} = ... # FIXME Constant\n"

    if ret and constants:
        ret += "\n"

    # Functions
    for name in sorted(functions):
        override = _check_override(prefix_name, name, local_overrides)
        if override:
            ret += override + "\n"
            continue

        ret += _build_function(
            current_namespace, name, functions[name], in_class, needed_namespaces
        )

    if ret and functions:
        ret += "\n"

    # Classes
    for name, obj in sorted(classes.items()):
        override = _check_override(prefix_name, name, local_overrides)
        if override:
            ret += override + "\n\n"
            continue

        full_name = _generate_full_name(prefix_name, name)

        class_ret = _gi_build_stub(
            repo,
            obj,
            current_namespace,
            _find_methods(obj),
            needed_namespaces,
            local_overrides,
            obj,
            full_name,
        )

        readable_props: list[GIRepository.BaseInfo] = []
        writable_props: list[GIRepository.BaseInfo] = []
        parents: list[str] = []
        fields: list[str] = []
        if hasattr(obj, "__info__"):
            object_info = obj.__info__  # type: ignore

            if isinstance(object_info, GI.StructInfo):
                for f in object_info.get_fields():
                    t = _type_to_python(
                        f.get_type_info(), current_namespace, needed_namespaces, True
                    )
                    n = f.get_name()
                    if n in dir(obj):
                        override = _check_override(full_name, n, local_overrides)
                        if override:
                            fields.append(override)
                        else:
                            fields.append(f"{n}: {t}")

            if isinstance(object_info, GI.ObjectInfo):
                p = object_info.get_parent()
                if p:
                    if current_namespace == p.get_namespace():
                        parents.append(f"{p.get_name()}")
                    else:
                        parents.append(f"{p.get_namespace()}.{p.get_name()}")
                        needed_namespaces.add(p.get_namespace())

                ifaces = object_info.get_interfaces()
                for i in ifaces:
                    if current_namespace == i.get_namespace():
                        parents.append(f"{i.get_name()}")
                    else:
                        parents.append(f"{i.get_namespace()}.{i.get_name()}")
                        needed_namespaces.add(i.get_namespace())

                for f in object_info.get_fields():
                    t = _type_to_python(
                        f.get_type_info(), current_namespace, needed_namespaces, True
                    )
                    n = f.get_name()
                    if n in dir(obj):
                        override = _check_override(full_name, n, local_overrides)
                        if override:
                            fields.append(override)
                        else:
                            fields.append(f"{n}: {t}")

                # Properties
                (rp, wp) = _object_get_props(repo, object_info)
                readable_props.extend(rp)
                writable_props.extend(wp)

            if isinstance(object_info, GI.InterfaceInfo):
                if current_namespace == "GObject":
                    parents.append("Object")
                else:
                    parents.append("GObject.GInterface")
                    needed_namespaces.add("GObject")

            if issubclass(obj, GObject.GBoxed):
                if current_namespace == "GObject":
                    parents.append("GBoxed")
                else:
                    parents.append("GObject.GBoxed")
                    needed_namespaces.add("GObject")

            if issubclass(obj, GObject.GPointer):
                if current_namespace == "GObject":
                    parents.append("GPointer")
                else:
                    parents.append("GObject.GPointer")
                    needed_namespaces.add("GObject")

        string_parents = ""
        if len(parents) > 0:
            string_parents = f"({', '.join(parents)})"

        if (
            not class_ret
            and len(fields) == 0
            and len(readable_props) == 0
            and len(writable_props) == 0
        ):
            ret += f"class {name}{string_parents}: ...\n"
        else:
            ret += f"class {name}{string_parents}:\n"

            # extracting docs
            doc = getattr(obj, "__doc__", "") or ""
            gdoc = getattr(obj, "__gdoc__", "") or ""

            txt = doc.strip() + "\n\n" + gdoc.strip()
            if not txt.isspace():
                txt = 'r"""\n' + txt.strip() + '\n"""' + "\n"
                txt = textwrap.indent(txt, "    ")
                ret += txt

        names: list[str]
        s: list[str]
        props_override = _check_override(full_name, "Props", local_overrides)
        if props_override:
            for line in props_override.splitlines():
                ret += "    " + line + "\n"
        elif len(readable_props) > 0 or len(writable_props) > 0:
            names = []
            s = []
            for p in itertools.chain(readable_props, writable_props):
                assert isinstance(p, GIRepository.PropertyInfo)
                local_name = p.get_name()
                assert local_name is not None
                n = local_name.replace("-", "_")
                if n in names:
                    # Avoid duplicates
                    continue
                names.append(n)
                ltype = p.get_type_info()
                t = _type_to_python(ltype, current_namespace, needed_namespaces, True)

                # Check getter/setter
                getter = p.get_getter()
                setter = p.get_setter()
                if getter and getter.may_return_null():
                    s.append(f"{n}: typing.Optional[{t}]")
                elif setter and not getter:
                    # If is writable only prop check if setter can accept NULL
                    arg_info = setter.get_arg(0)
                    if arg_info.may_be_null():
                        s.append(f"{n}: typing.Optional[{t}]")
                    else:
                        s.append(f"{n}: {t}")
                else:
                    s.append(f"{n}: {t}")

            separator = "\n        "
            ret += f"    class Props:\n        {separator.join(s)}\n"
            ret += "    props: Props = ...\n"

        for field in fields:
            ret += f"    {field} = ...\n"

        class_constructor_override = _check_override(
            full_name, "__init__", local_overrides
        )
        if class_constructor_override:
            for line in class_constructor_override.splitlines():
                ret += "    " + line + "\n"
        elif len(writable_props) > 0:
            names = []
            s = []
            for p in writable_props:
                assert isinstance(p, GIRepository.PropertyInfo)
                local_name = p.get_name()
                assert local_name is not None
                n = local_name.replace("-", "_")
                if n in names:
                    # Avoid duplicates
                    continue
                names.append(n)
                ltype = p.get_type_info()
                t = _type_to_python(ltype, current_namespace, needed_namespaces)
                setter = p.get_setter()
                if setter:
                    arg_info = setter.get_arg(0)
                    if arg_info.may_be_null():
                        s.append(f"{n}: typing.Optional[{t}] = ...")
                    else:
                        s.append(f"{n}: {t} = ...")
                else:
                    s.append(f"{n}: {t} = ...")

            separator = ",\n                 "
            ret += f"    def __init__(self, {separator.join(s)}) -> None: ...\n"

        for line in class_ret.splitlines():
            ret += "    " + line + "\n"

        ret += "\n"

    # Flags
    for name, obj in sorted(flags.items()):
        override = _check_override(prefix_name, name, local_overrides)
        if override:
            ret += override + "\n\n"
            continue

        full_name = _generate_full_name(prefix_name, name)

        if current_namespace == "GObject":
            if name != "GFlags":
                base = "GFlags"
            else:
                base = ""
        else:
            needed_namespaces.add("GObject")
            base = "GObject.GFlags"

        ret += f"class {name}({base}):\n"
        for key in sorted(vars(obj)):
            if key.startswith("__") or key[0].isdigit():
                continue

            override = _check_override(full_name, key, local_overrides)
            if override:
                for line in override.splitlines():
                    ret += "    " + line + "\n"
                continue

            o = getattr(obj, key)
            if isinstance(o, GI.FunctionInfo):
                function_ret = _build_function(
                    current_namespace, key, o, obj, needed_namespaces
                )
                for line in function_ret.splitlines():
                    ret += "    " + line + "\n"
            else:
                ret += f"    {key}: {o.__class__.__name__} = ...\n"
        ret += "\n"

    # Enums
    for name, obj in sorted(enums.items()):
        override = _check_override(prefix_name, name, local_overrides)
        if override:
            ret += override + "\n\n"
            continue

        full_name = _generate_full_name(prefix_name, name)

        if current_namespace == "GObject":
            if name != "GEnum":
                base = "GEnum"
            else:
                base = ""
        else:
            needed_namespaces.add("GObject")
            base = "GObject.GEnum"

        # Some Enums can be empty in the end
        ret += f"class {name}({base}):\n"
        length_before = len(ret)
        for key in sorted(vars(obj)):
            if key.startswith("__") or key[0].isdigit():
                continue

            override = _check_override(full_name, key, local_overrides)
            if override:
                for line in override.splitlines():
                    ret += "    " + line + "\n"
                continue

            o = getattr(obj, key)
            if isinstance(o, GI.FunctionInfo):
                function_ret = _build_function(
                    current_namespace, key, o, obj, needed_namespaces
                )
                for line in function_ret.splitlines():
                    ret += "    " + line + "\n"
            else:
                ret += f"    {key}: {o.__class__.__name__} = ...\n"

        if len(ret) == length_before:
            # No attributes where found
            ret += "    ...\n"

        ret += "\n"

    return ret


def _find_methods(obj: Type[Any]) -> list[str]:
    mro = inspect.getmro(obj)
    main_name = _get_gname(mro[0])

    all_attrs = set(dir(obj))
    other_attrs: set[str] = set()
    for o in mro[1:]:
        name = _get_gname(o)
        if name == main_name:
            continue
        other_attrs.update(dir(o))

    obj_attrs = all_attrs - other_attrs

    # Search for overridden methods
    if hasattr(obj, "__info__"):
        obj_info = obj.__info__  # type: ignore
        if isinstance(obj_info, (GI.ObjectInfo, GI.StructInfo)):
            methods = obj_info.get_methods()
            for m in methods:
                name = m.get_name()
                if name in dir(obj) and name not in obj_attrs:
                    obj_attrs.add(name)

    return sorted(list(obj_attrs))


def _get_gname(obj: Type[Any]) -> Optional[str]:
    if not hasattr(obj, "__gtype__"):
        return None
    return obj.__gtype__.name  # type: ignore


def start(module: str, version: str, local_overrides: dict[str, str]) -> str:
    """Start."""
    gi.require_version(module, version)
    m = importlib.import_module(f".{module}", "gi.repository")
    repo = GIRepository.Repository.new()
    repo.require(module, version, cast(GIRepository.RepositoryLoadFlags, 0))
    return _build(repo, m, module, local_overrides)


def main_func() -> None:
    """Execute main."""
    parser = argparse.ArgumentParser(description="Generate module stubs")
    parser.add_argument("module", type=str, help="Gdk, Gtk, ...")
    parser.add_argument("version", type=str, help="3.0, 4.0, ...")
    parser.add_argument(
        "-u", dest="update", type=str, help="Stub file to update e.g. -u Gdk.pyi"
    )

    args = parser.parse_args()

    if args.update:
        overrides: dict[str, str] = {}
        try:
            with open(args.update, "r", encoding="utf-8") as file:
                overrides = parse_pyi.parse(file.read())
        except FileNotFoundError:
            pass
        output = start(args.module, args.version, overrides)
        print("Running with this overrides:")
        pprint.pprint(overrides)
        with open(args.update, "w+", encoding="utf-8") as file:
            file.write(output)
    else:
        print(start(args.module, args.version, {}))


if __name__ == "__main__":
    main_func()
