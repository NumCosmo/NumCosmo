from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import powerlaw_free_emissivity,rename_machine,rename_replace,norename

cdef extern c_parametric *gal545_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_free_emissivity_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_free_emissivity_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class powerlaw_free_emissivity_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_EE_init

cdef class powerlaw_free_emissivity_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_TE_init

def galf_P_rename_func(v,rups):
  if v.startswith("galf"):
    rv = v.replace("galf","pwfe")
    rups[v]=rv

galf_TE = rename_machine(powerlaw_free_emissivity_TE,{"galf_TE_l2_norm":"1","galf_TE_l_pivot":"500","galf_TE_index":"-2.4"},galf_P_rename_func)
galf_EE = rename_machine(powerlaw_free_emissivity_EE,{"galf_EE_l2_norm":"1","galf_EE_l_pivot":"500","galf_EE_index":"-2.4"},galf_P_rename_func)

def ps_P_rename_func(v,rups):
  if v.startswith("ps"):
    rv = v.replace("ps","pwfe")
    rups[v]=rv

ps_TE = rename_machine(powerlaw_free_emissivity_TE,{"ps_TE_l2_norm":"1","ps_TE_l_pivot":"3000","ps_TE_index":"0"},ps_P_rename_func)
ps_EE = rename_machine(powerlaw_free_emissivity_EE,{"ps_EE_l2_norm":"1","ps_EE_l_pivot":"3000","ps_EE_index":"0"},ps_P_rename_func)


cdef class gal545(parametric):
  def __cinit__(self):
    self.initfunc = <void*> gal545_init;


component_list = ["gal545","powerlaw_free_emissivity_EE","powerlaw_free_emissivity_TE","galf_TE","galf_EE","ps_TE","ps_EE"]


cdef extern c_parametric *pointsource_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class pointsource(parametric):
  def __cinit__(self):
    self.initfunc = <void*>pointsource_init

component_list += ["pointsource"]


cdef extern double c_sz_spectrum "sz_spectrum" (double nu, double nu0)
cdef extern c_parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ksz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric  *gibXsz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
cdef extern c_parametric *gcib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

def sz_spectrum(nu,nu0=143.0):
  return c_sz_spectrum(<double>nu,<double>nu0)

pname = "rel2015"

cdef class sz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> sz_init;
    self.template_name = "tsz_143_eps0.50.dat[1]"
    self.plugin_name = pname

cdef class ksz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ksz_init;
    self.template_name = "ksz_fromcamspec.dat"
    self.plugin_name = pname

cdef class gcib(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> gcib_init;
    self.template_name = "cib_1h_2h_100_353_Jsr-1_PS_2014_09.dat"
    self.plugin_name = pname

cdef class gibXsz(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> gibXsz_init;
    self.template_name = ["sz_x_cib_template.dat[1]"]
    self.plugin_name = pname
    

cib_1h_2h_sept14 = rename_machine(gcib,{},rename_replace("gib","cib"))
sz_color = rename_machine(sz,{"sz_color_143_to_143":"0.975","sz_color_100_to_143":"0.981"},norename)

component_list += ["sz_color","gibXsz","gcib","sz","ksz","cib_1h_2h_sept14"]

cdef extern c_parametric *bleak_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cnoise_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

cdef class bleak(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> bleak_init;
    self.template_name = "sky_template_v15_F100_143_217_353.dat"
    self.plugin_name = pname

cdef class cnoise(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> cnoise_init;
    self.template_name = "cnoise_F100_143_217_353_v17.dat"
    self.plugin_name = pname

cnoise_v17 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_v17.dat")
bleak_v15 = rename_machine(bleak,{},norename,data_file="sky_template_v15_F100_143_217_353.dat")    

component_list += ["bleak","cnoise","cnoise_v17","bleak_v15"]
 

