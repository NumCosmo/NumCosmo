-- Nouns used to replace a crossref reference once its target is gone.
local crossref_nouns = {
  eq = "equation",
  fig = "figure",
  lst = "listing",
  sec = "section",
  tbl = "table",
  thm = "theorem",
}

-- Flattening a FloatRefTarget below drops its crossref anchor, so a surviving
-- "@tbl-foo" cannot resolve: Quarto warns and the notebook renders a literal
-- "?@tbl-foo". Replace such references with plain text.
local function replace_crossref(cite)
  if #cite.citations ~= 1 then
    return nil
  end

  local prefix = cite.citations[1].id:match("^(%a+)%-")
  local noun = prefix and crossref_nouns[prefix]

  if not noun then
    return nil
  end

  return { pandoc.Str("the"), pandoc.Space(), pandoc.Str(noun) }
end

function Pandoc(doc)
  if not quarto.doc.isFormat("ipynb") then
    return doc
  end

  for i, b in ipairs(doc.blocks) do
    if b.t == "Div" then
      local a = b.attr and b.attr.attributes or {}

      if a["__quarto_custom_type"] == "FloatRefTarget" then
        local code = nil

        b:walk {
          CodeBlock = function(cb)
            if not code then
              code = cb
            end
          end
        }

        if code then
          quarto.log.output("CLONING EXEC CELL")

          doc.blocks[i] = pandoc.Div(
            { code },
            pandoc.Attr(
              "",
              { "cell" },
              {
                execution_count = "0"
              }
            )
          )
        end
      end
    end
  end

  return doc:walk { Cite = replace_crossref }
end
