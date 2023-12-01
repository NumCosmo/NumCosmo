#!/usr/bin/env python3
#
# addmathjax.py
#
# Wed Nov 29 10:52:10 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# addmathjax.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Add MathJax script tag to HTML files in-place."""

import sys
import os
from bs4 import BeautifulSoup


def add_mathjax_script_inplace(filename):
    """Add MathJax script tag to HTML files in-place."""
    with open(filename, "r", encoding="utf-8") as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, "html.parser")
    head_tag = soup.head

    link_tag = soup.new_tag(
        "link",
        rel="stylesheet",
        href="container.css",
        type="text/css",
    )
    head_tag.append(link_tag)

    script_conf_tag = soup.new_tag(
        "script",
        type="text/x-mathjax-config",
    )
    script_conf_tag.append(
        """
//<![CDATA[
MathJax.Hub.Config({"HTML-CSS": { preferredFont: "TeX", availableFonts: ["STIX","TeX"],
linebreaks: { automatic:true }, EqnChunk: (MathJax.Hub.Browser.isMobile ? 10 : 50) },
tex2jax: { inlineMath: [ ["$", "$"], ["\\\\(","\\\\)"] ], displayMath: [ ["$$","$$"], 
["\\[", "\\]"] ], processEscapes: true, ignoreClass: "tex2jax_ignore|dno" },
TeX: {  noUndefined: { 
    attributes: { mathcolor: "red", mathbackground: "#FFEEEE", mathsize: "90%" 
    } },
equationNumbers: { autoNumber: "AMS" } },
messageStyle: "none"
});
//]]>
        """
    )
    head_tag.append(script_conf_tag)

    script_tag = soup.new_tag(
        "script",
        type="text/javascript",
        src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/"
        "MathJax.js?config=TeX-AMS_HTML",
    )
    head_tag.append(script_tag)

    modified_html = str(soup)

    with open(filename, "w", encoding="utf-8") as file:
        file.write(modified_html)


def run_addmathjax():
    """Run add in-place MathJax script tag to HTML files in a directory."""

    if len(sys.argv) != 2:
        print("Usage: python addmathjax.py <html_directory>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: directory '{directory}' does not exist.")
        sys.exit(1)

    for filename in os.listdir(directory):
        if filename.endswith(".html"):
            filename0 = os.path.join(directory, filename)
            add_mathjax_script_inplace(filename0)


if __name__ == "__main__":
    run_addmathjax()
