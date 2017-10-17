#!/usr/bin/perl

use Glib::Object::Introspection;

Glib::Object::Introspection->setup (
  basename => 'NumCosmoMath',
  version => '1.0',
  package => 'Ncm');

Glib::Object::Introspection->setup (
  basename => 'NumCosmo',
  version => '1.0',
  package => 'Nc');

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm::cfg_init ();

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#
my $cosmo = Nc::HICosmo->new_from_name ("Nc::HICosmo", "NcHICosmoDEXcdm");

#
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
$dist = new Nc::Distance (2.0);

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remember to use the _orig_ version to set the original
#  parameters when a reparametrization is used.
#
my $h0 = Glib::Object::Introspection->convert_sv_to_enum ("Nc::HICosmoDESParams", "h0");


#
# C-like
#
$cosmo->orig_param_set ($h0,        70.0);

print $cosmo->set_property ("H0", 78.0);

#
# OO-like
#
$cosmo{H0}      = 72.0;
$cosmo{Omegab}  = 0.05;
$cosmo{Omegac}  = 0.25;
$cosmo{Omegax}  = 0.70;
$cosmo{Tgamma0} = 2.72;
$cosmo{ns}      = 1.0;
$cosmo{sigma8}  = 0.9;
$cosmo{w}       = -1.0;

#
#  Printing the parameters used.
#
print "# Model parameters: ";
$cosmo->params_log_all ();

#
#  Printing some distances up to redshift 1.0.
#
for (my $i = 0; $i < 10; $i++) {
  my $z = 1.0 / 9.0 * $i;
  $cd = Ncm::C::hubble_radius () * $dist->comoving ($cosmo, $z);
  print "$z $cd\n";

}

