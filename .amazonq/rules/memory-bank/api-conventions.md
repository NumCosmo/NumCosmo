# NumCosmo API Conventions

## Naming Conventions

### Namespace Prefixes
- **NCM** - NumCosmo Math (general mathematical utilities)
- **NC** - NumCosmo (cosmology-specific implementations)

### Class Naming Pattern
- Format: `[prefix]_[module]_[class]`
- Examples: `NcmModel`, `NcHICosmo`, `NcDistance`, `NcmSpline`

### Function Naming Pattern
- Format: `[class_name]_[method_name]`
- Examples: `ncm_model_new()`, `nc_hicosmo_E()`, `nc_distance_comoving()`
- Constructor: `_new()` or `_new_from_*()`
- Destructor: `_free()` or use `g_object_unref()`
- Getter: `_get_*()` or `_peek_*()`
  - `_get_*()` increases reference count (caller must unref)
  - `_peek_*()` returns borrowed reference without increasing count (use when object won't be saved)
- Setter: `_set_*()`

### File Naming
- Header files: `[class_name].h`
- Implementation files: `[class_name].c`
- Example: `nc_distance.h`, `nc_distance.c`

## GObject Patterns

### Object Lifecycle
- Objects created with `_new()` functions
- Reference counting via `g_object_ref()` / `g_object_unref()`
- Use `_clear()` macros for safe unreferencing: `ncm_model_clear(&model)`
- Use `_free()` for non-GObject structures

### Properties
- Access via `g_object_get()` / `g_object_set()`
- Property names use kebab-case: "property-name"
- Can also use dedicated getter/setter functions

### Type System
- Type registration: `NC_TYPE_*` or `NCM_TYPE_*`
- Type checking: `NC_IS_*()` or `NCM_IS_*()`
- Type casting: `NC_*()` or `NCM_*()`

## Memory Management

### Ownership Rules
- Functions returning `_new()` transfer full ownership
- Functions with `_get_` prefix return owned references (increment ref count)
- Functions with `_peek_` prefix return borrowed references (no ref count change)
  - Use `_peek_*()` in C when temporarily using an object without saving it
  - Caller must not unref a peeked reference
- Functions with `_dup_` prefix return owned copies
- Arrays and strings follow GLib conventions

### Reference Counting
- Always unref objects when done
- Use `_clear()` macros to avoid double-free
- Parent objects hold references to children

## Parameter Conventions

### Input Parameters
- Const pointers for read-only objects: `const NcmModel *model`
- Non-const for objects being modified: `NcmModel *model`
- Output parameters use double pointers: `gdouble *result`

### Return Values
- Success/failure: boolean `gboolean`
- Floating point: `gdouble` (not `double`)
- Integers: `gint`, `guint` (not `int`)
- Sizes: `gsize` or `guint`

## Error Handling

### GError Pattern
```c
gboolean function_name (NcmObject *obj, GError **error);
```
- Return FALSE on error, TRUE on success
- Set error details via `g_set_error()`
- Check for NULL error pointer before setting

### Assertions
- Use `g_assert()` for internal consistency checks
- Use `g_return_if_fail()` for parameter validation
- Use `g_return_val_if_fail()` when returning a value

## Documentation

### GTK-Doc Format
- Brief description on first line
- Detailed description in following paragraphs
- Parameter documentation: `@param_name: description`
- Return value: `Returns: description`
- Since version: `Since: X.Y.Z`

### Example
```c
/**
 * nc_distance_comoving:
 * @dist: a #NcDistance
 * @z: redshift
 *
 * Computes the comoving distance to redshift @z.
 *
 * Returns: comoving distance in Mpc
 */
```

## Python Bindings

### Automatic Introspection
- GObject Introspection generates Python bindings automatically
- Python naming converts underscores to camelCase for methods
- C: `nc_distance_comoving()` → Python: `dist.comoving()`

### Python Conventions
- Import: `from numcosmo_py import Nc, Ncm`
- Object creation: `dist = Nc.Distance.new(z_max)`
- Method calls: `result = dist.comoving(z)`

## Thread Safety

### General Rules
- Most objects are NOT thread-safe by default
- Use separate instances per thread
- Some objects support explicit thread-safe operations
- Document thread-safety in class description

### Parallel Computing
- MPI support for distributed computing
- OpenMP for shared memory parallelization
- Use appropriate locks when sharing objects
