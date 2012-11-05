#!/usr/bin/perl -w


use strict;

use XML::LibXML;

my $file = $ARGV[0];
my $nsURI = "http://docbook.org/ns/docbook";
my $nsp = "docbooka";

my $parser = XML::LibXML->new ();
my $doc = $parser->load_xml (location => $file);
$doc->removeInternalSubset ();
my $dtd = $doc->createInternalSubset (
  "article", 
  "-//OASIS//DTD DocBook XML V4.5//EN", 
  "http://www.oasis-open.org/docbook/xml/4.5/docbookx.dtd") or die $!;

my $root = $doc->documentElement();
my $newroot = $doc->createElement("chapter");
#$newroot->setNamespace($nsURI, undef, 0);

my $xpc = XML::LibXML::XPathContext->new ($doc);
$xpc->registerNs ($nsp, $nsURI);

foreach my $eq ($xpc->findnodes("//$nsp:equation")) {
  foreach my $note ($xpc->findnodes("following-sibling::$nsp:note[\@role='equation']", $eq)) {
    $note->parentNode->removeChild ($note);
  }
}

foreach my $note_eq ($xpc->findnodes("//$nsp:note[\@role='equation']")) {
  $note_eq->parentNode->removeChild ($note_eq);
}

foreach my $link ($xpc->findnodes("//$nsp:link")) {
  if ($link->hasAttributeNS ("http://www.w3.org/1999/xlink", "href")) {
    my $href = $link->getAttributeNS ("http://www.w3.org/1999/xlink", "href");
    $link->removeAttributeNS ("http://www.w3.org/1999/xlink", "href");
    $href =~ s/^\#//;
    $link->setAttribute("linkend", $href);
  }
  print $link->getAttributeNodeNS ("http://www.w3.org/1999/xlink", "href"), "\n";
  print $link->toString;
}

my @secs = $xpc->findnodes("/$nsp:article/$nsp:section");

my @titles = $xpc->findnodes("//$nsp:note[\@role=\"title\"]");

if (@titles == 1) {
  my $title = $doc->createElement("title");
  $title->appendText($titles[0]->textContent);
  $newroot->appendChild ($title);
}

foreach my $sec (@secs) {
  my $newsec = $parser->load_xml(string => $sec->toString);
  
  foreach my $abbrev ($newsec->findnodes ("//title/abbrev")) {
    $abbrev->parentNode->removeChild ($abbrev);
  }

  $newroot->appendChild ($newsec->documentElement);
}

$doc->setDocumentElement($newroot);

open ARQ, ">$file";
print ARQ $doc->toString (1 );
close ARQ;
#print $doc->toString (1 );
