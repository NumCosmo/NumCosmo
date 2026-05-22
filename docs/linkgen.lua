-- Auto-link NumCosmo API symbols in Quarto text.
-- Supported forms:
--   [[numcosmo|NcHIPertAdiab]]
--   [[numcosmo-math|NcmCSQ1D]]
--   [[ncm_csq1d_eval_delta_theta_at]]

local namespaces = {
    ["numcosmo"] = {
        href_base = "reference/numcosmo",
        rel_ref = "numcosmo",
    },
    ["numcosmo-math"] = {
        href_base = "reference/numcosmo-math",
        rel_ref = "numcosmo-math",
    },
    ["math"] = {
        href_base = "reference/numcosmo-math",
        rel_ref = "numcosmo-math",
    },
}

local type_prefix = {
    alias = "alias",
    bitfield = "flags",
    callback = "callback",
    class = "class",
    class_method = "class_method",
    constant = "const",
    ctor = "ctor",
    domain = "error",
    enum = "enum",
    ["function"] = "func",
    ["function_macro"] = "func",
    interface = "iface",
    method = "method",
    property = "property",
    record = "struct",
    type_func = "type_func",
    vfunc = "vfunc",
}

local cache = {}

local function path_join(a, b)
    if a == "" then
        return b
    end
    if a:sub(-1) == "/" then
        return a .. b
    end
    return a .. "/" .. b
end

local function dirname(path)
    if path == nil or path == "" then
        return ""
    end
    local d = path:match("^(.*)/[^/]+$")
    if d == nil then
        return ""
    end
    return d
end

local function push_unique(tbl, value)
    for _, v in ipairs(tbl) do
        if v == value then
            return
        end
    end
    table.insert(tbl, value)
end

local function shell_quote(s)
    return "'" .. tostring(s):gsub("'", "'\\''") .. "'"
end

local function list_direct_subdirs(path)
    local dirs = {}
    local cmd = "find " .. shell_quote(path) .. " -mindepth 1 -maxdepth 1 -type d -print 2>/dev/null"
    local p = io.popen(cmd)
    if p == nil then
        return dirs
    end
    for line in p:lines() do
        if line ~= "" then
            table.insert(dirs, line)
        end
    end
    p:close()
    return dirs
end

local function script_dir()
    local src = debug.getinfo(1, "S").source
    if src ~= nil and src:sub(1, 1) == "@" then
        return dirname(src:sub(2))
    end
    return ""
end

local function project_roots()
    local roots = { "." }

    local sdir = script_dir()
    if sdir ~= "" then
        push_unique(roots, sdir)
        local parent = dirname(sdir)
        if parent ~= "" then
            push_unique(roots, parent)
        end
    end

    local input = ""
    if PANDOC_STATE ~= nil and PANDOC_STATE.input_files ~= nil and PANDOC_STATE.input_files[1] ~= nil then
        input = PANDOC_STATE.input_files[1]
    end

    if input ~= "" then
        local docs_root = input:match("^(.*)/docs/.*$")
        if docs_root ~= nil and docs_root ~= "" then
            push_unique(roots, docs_root)
            push_unique(roots, path_join(docs_root, "docs"))
        end
    end

    return roots
end

local function lookup_candidates(ns)
    local rel = ns.rel_ref
    local out = {}

    for _, root in ipairs(project_roots()) do
        push_unique(out, path_join(root, path_join("reference", rel)))
        push_unique(out, path_join(root, path_join("docs/reference", rel)))

        for _, child in ipairs(list_direct_subdirs(root)) do
            push_unique(out, path_join(child, path_join("docs/reference", rel)))
        end
    end

    return out
end

local function file_exists(path)
    local f = io.open(path, "r")
    if f == nil then
        return false
    end
    f:close()
    return true
end

local function read_all(path)
    local f = io.open(path, "r")
    if f == nil then
        return nil
    end
    local content = f:read("*a")
    f:close()
    return content
end

local function add_unique(tbl, key, value)
    local cur = tbl[key]

    if cur == nil then
        tbl[key] = value
        return
    end

    if cur == value then
        return
    end

    -- Mark ambiguous names so we do not create a potentially wrong link.
    tbl[key] = false
end

local function symbol_to_href(base, sym)
    local prefix = type_prefix[sym.type]
    if prefix == nil then
        return nil
    end

    if sym.type == "content" then
        return nil
    end

    if sym.type == "ctor" or sym.type == "method" or sym.type == "property" or sym.type == "type_func" or sym.type == "vfunc" then
        if sym.type_name == nil or sym.name == nil then
            return nil
        end
        return string.format("%s/%s.%s.%s.html", base, prefix, sym.type_name, sym.name)
    end

    if sym.type == "class_method" then
        local owner = sym.struct_for or sym.type_name
        if owner == nil or sym.name == nil then
            return nil
        end
        return string.format("%s/%s.%s.%s.html", base, prefix, owner, sym.name)
    end

    if sym.name == nil then
        return nil
    end

    return string.format("%s/%s.%s.html", base, prefix, sym.name)
end

local function index_path(base)
    return base .. "/index.json"
end

local function pick_existing_base(candidates)
    for _, base in ipairs(candidates) do
        if file_exists(index_path(base)) then
            return base
        end
    end

    -- Fallback to first candidate to keep deterministic href generation.
    return candidates[1]
end

local function load_namespace(key)
    local ns = namespaces[key]
    if ns == nil then
        return nil
    end

    local lookup_base = pick_existing_base(lookup_candidates(ns))
    local href_base = ns.href_base

    if cache[key] ~= nil then
        return cache[key]
    end

    local idx = {
        by_ident = {},
        by_ctype = {},
        by_name = {},
    }

    local p = index_path(lookup_base)
    local raw = read_all(p)
    if raw == nil then
        cache[key] = idx
        return idx
    end

    local ok, doc = pcall(pandoc.json.decode, raw)
    if not ok or type(doc) ~= "table" or type(doc.symbols) ~= "table" then
        cache[key] = idx
        return idx
    end

    for _, sym in ipairs(doc.symbols) do
        local href = symbol_to_href(href_base, sym)
        if href ~= nil then
            if sym.ident ~= nil then
                add_unique(idx.by_ident, sym.ident, href)
            end
            if sym.ctype ~= nil then
                add_unique(idx.by_ctype, sym.ctype, href)
            end
            if sym.name ~= nil then
                add_unique(idx.by_name, sym.name, href)
            end
        end
    end

    cache[key] = idx
    return idx
end

local function normalize_target(target)
    if target == nil then
        return nil
    end

    local lowered = target:lower()
    if lowered == "numcosmo" then
        return "numcosmo"
    end
    if lowered == "numcosmo-math" or lowered == "math" or lowered == "ncm" then
        return "numcosmo-math"
    end

    return nil
end

local function resolve_in_ns(ns, symbol)
    local idx = load_namespace(ns)
    if idx == nil then
        return nil
    end

    local by_ident = idx.by_ident[symbol]
    if type(by_ident) == "string" then
        return by_ident
    end

    local by_ctype = idx.by_ctype[symbol]
    if type(by_ctype) == "string" then
        return by_ctype
    end

    local by_name = idx.by_name[symbol]
    if type(by_name) == "string" then
        return by_name
    end

    -- Allow type-like shortcuts, e.g. NcHIPertAdiab -> HIPertAdiab.
    local short = symbol:gsub("^Nc", ""):gsub("^Ncm", "")
    if short ~= symbol then
        local short_name = idx.by_name[short]
        if type(short_name) == "string" then
            return short_name
        end
    end

    return nil
end

local function resolve_symbol(target, symbol)
    local ns = normalize_target(target)

    if ns ~= nil then
        return resolve_in_ns(ns, symbol)
    end

    local first = resolve_in_ns("numcosmo", symbol)
    local second = resolve_in_ns("numcosmo-math", symbol)

    if first ~= nil and second == nil then
        return first
    end
    if second ~= nil and first == nil then
        return second
    end

    return nil
end

local function root_prefix()
    local input = ""
    if PANDOC_STATE ~= nil and PANDOC_STATE.input_files ~= nil and PANDOC_STATE.input_files[1] ~= nil then
        input = PANDOC_STATE.input_files[1]
    end

    if input == "" then
        return ""
    end

    -- Prefer the path relative to the docs root when an absolute input path is provided.
    local rel = input
    local after_docs = input:match(".*/docs/(.*)")
    if after_docs ~= nil and after_docs ~= "" then
        rel = after_docs
    end

    local dir = rel:match("^(.*)/[^/]+$")
    if dir == nil or dir == "" then
        return ""
    end

    local depth = 0
    for _ in dir:gmatch("[^/]+") do
        depth = depth + 1
    end

    return string.rep("../", depth)
end

local function build_link(target, symbol)
    local href = resolve_symbol(target, symbol)
    if href == nil then
        return nil
    end

    local url = root_prefix() .. href
    return pandoc.Link(pandoc.Code(symbol), url)
end

local function parse_token(text)
    local pfx, target, symbol, sfx = text:match("(.-)%[%[(.-)|(.-)%]%](.*)")
    if target ~= nil and symbol ~= nil then
        return pfx, target, symbol, sfx
    end

    local pfx2, symbol2, sfx2 = text:match("(.-)%[%[(.-)%]%](.*)")
    if symbol2 ~= nil then
        return pfx2, nil, symbol2, sfx2
    end

    return nil
end

return {
    {
        Inline = function(el)
            if el.t ~= "Str" then
                return el
            end

            local prefix, target, symbol, suffix = parse_token(el.text)
            if symbol == nil then
                return el
            end

            local link = build_link(target, symbol)
            if link == nil then
                return el
            end

            local out = {}
            if prefix ~= "" then
                table.insert(out, pandoc.Str(prefix))
            end

            table.insert(out, link)

            if suffix ~= "" then
                table.insert(out, pandoc.Str(suffix))
            end

            return out
        end,
    },
}
