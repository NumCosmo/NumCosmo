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

  return doc
end
