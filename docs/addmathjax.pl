#!/usr/bin/perl -w

use strict;

use XML::LibXML;

my $parser = XML::LibXML->new (recover => 1);

while (my $file = shift(@ARGV)) {

  my $doc = $parser->load_html (location => $file);

  my $dtd = $doc->createExternalSubset ("html", "-//W3C//DTD XHTML 1.0 Transitional//EN", "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd");
  
  $doc->setInternalSubset($dtd);

  my $root = $doc->documentElement();

  my @heads = $root->getChildrenByTagName("head");
  
  my $mathjax_chunk = <<'ENDSCRIPT';
<link rel="stylesheet" href="container.css" type="text/css"></link>
<script type="text/x-mathjax-config">
//<![CDATA[
                MathJax.Hub.Config({"HTML-CSS": { preferredFont: "TeX", availableFonts: ["STIX","TeX"], linebreaks: { automatic:true }, EqnChunk: (MathJax.Hub.Browser.isMobile ? 10 : 50) },
                                    tex2jax: { inlineMath: [ ["$", "$"], ["\\\\(","\\\\)"] ], displayMath: [ ["$$","$$"], ["\\[", "\\]"] ], processEscapes: true, ignoreClass: "tex2jax_ignore|dno" },
                                    TeX: {  noUndefined: { attributes: { mathcolor: "red", mathbackground: "#FFEEEE", mathsize: "90%" } },
                                            equationNumbers: { autoNumber: "AMS" } },
                                    messageStyle: "none"
                });
//]]>
</script>
<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML"></script> 
ENDSCRIPT
  my $mathjax = $parser->parse_balanced_chunk ($mathjax_chunk);
  $heads[0]->appendChild ($mathjax);

  my @bodys = $root->getChildrenByTagName("body");
  my @childnodes = $bodys[0]->childNodes();

  my $div1 = XML::LibXML::Element->new ("div");
  $div1->setAttribute ('id', 'container');
  $div1->setAttribute ('class', 'gtk-doc');
  my $div2 = XML::LibXML::Element->new ("div");
  $div2->setAttribute ('class', 'content');
  $div1->appendChild ($div2);
  
  foreach my $node (@childnodes) {
    $div2->appendChild ($node);
    $bodys[0]->removeChild ($node);
  }
  
  $bodys[0]->appendChild ($div1);

  open ARQ, ">$file";
  print ARQ $doc->toString (2);
  close ARQ;
}


