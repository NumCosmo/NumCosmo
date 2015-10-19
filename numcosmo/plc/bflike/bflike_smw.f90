! Code for evaluation of CMB likelihood in pixel space
! Author: L. Colombo <colombo@usc.edu> 
! Original code based on Tegmark & de Oliveira-Costa  arXiv:astro-ph/0012120
! This version implements a low rank update to speed up computations
!
!**************************************************************************
!* THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR   *    
!* IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED         *    
!* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE *    
!* DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,    *    
!* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES     *    
!* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR     *    
!* SERVICES; LOSS OF USE, DATA, OR PROFITS; BUSINESS INTERRUPTION; OR     *    
!* UNFOUNDED CONCLUSIONS ON THE NATURE OF LIFE, THE UNIVERSE AND          *
!* EVERYTHING) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN  * 
!* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)*
!* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF *
!* THE POSSIBILITY OF SUCH DAMAGE.                                        *    
!**************************************************************************

module fiducial_cls_mod_smw

  implicit none

!container for a fiducial set of cls, used to fill the covariance if not enough cls are provided
  private

  type fiducial_cls

     real(kind=8) :: cls(2:130,6)

   contains
     
     procedure :: init

  end type fiducial_cls

  type(fiducial_cls) ,save :: fid

  public fid

contains

  subroutine init(self)
    class(fiducial_cls) ,intent(inout) :: self

    self%cls(002,:) = [  1.095000d+03,  3.404900d-02,  1.713500d-06,  2.637800d+00, 0.d0, 0.d0]
    self%cls(003,:) = [  1.024300d+03,  4.671000d-02,  3.421400d-06,  3.085800d+00, 0.d0, 0.d0]
    self%cls(004,:) = [  9.598900d+02,  4.349400d-02,  5.689800d-06,  2.977700d+00, 0.d0, 0.d0]
    self%cls(005,:) = [  9.117200d+02,  3.135300d-02,  8.511700d-06,  2.594600d+00, 0.d0, 0.d0]
    self%cls(006,:) = [  8.780600d+02,  1.870400d-02,  1.187800d-05,  2.143800d+00, 0.d0, 0.d0]
    self%cls(007,:) = [  8.559000d+02,  1.010800d-02,  1.577800d-05,  1.720000d+00, 0.d0, 0.d0]
    self%cls(008,:) = [  8.417000d+02,  5.913000d-03,  2.020100d-05,  1.361900d+00, 0.d0, 0.d0]
    self%cls(009,:) = [  8.337700d+02,  4.357500d-03,  2.513400d-05,  1.100100d+00, 0.d0, 0.d0]
    self%cls(010,:) = [  8.305000d+02,  3.705500d-03,  3.056300d-05,  9.396400d-01, 0.d0, 0.d0]
    self%cls(011,:) = [  8.307300d+02,  3.200100d-03,  3.647400d-05,  8.678200d-01, 0.d0, 0.d0]
    self%cls(012,:) = [  8.335100d+02,  2.835000d-03,  4.285100d-05,  8.591900d-01, 0.d0, 0.d0]
    self%cls(013,:) = [  8.380700d+02,  2.707400d-03,  4.968100d-05,  8.883400d-01, 0.d0, 0.d0]
    self%cls(014,:) = [  8.442200d+02,  2.778700d-03,  5.694700d-05,  9.403700d-01, 0.d0, 0.d0]
    self%cls(015,:) = [  8.519200d+02,  2.976400d-03,  6.463600d-05,  1.003000d+00, 0.d0, 0.d0]
    self%cls(016,:) = [  8.611100d+02,  3.243600d-03,  7.273300d-05,  1.066000d+00, 0.d0, 0.d0]
    self%cls(017,:) = [  8.715800d+02,  3.586400d-03,  8.122500d-05,  1.127800d+00, 0.d0, 0.d0]
    self%cls(018,:) = [  8.831000d+02,  4.026700d-03,  9.009900d-05,  1.189000d+00, 0.d0, 0.d0]
    self%cls(019,:) = [  8.954500d+02,  4.586100d-03,  9.934500d-05,  1.250100d+00, 0.d0, 0.d0]
    self%cls(020,:) = [  9.084100d+02,  5.286600d-03,  1.089500d-04,  1.311700d+00, 0.d0, 0.d0]
    self%cls(021,:) = [  9.217700d+02,  6.146500d-03,  1.189200d-04,  1.374100d+00, 0.d0, 0.d0]
    self%cls(022,:) = [  9.355000d+02,  7.171600d-03,  1.292300d-04,  1.436700d+00, 0.d0, 0.d0]
    self%cls(023,:) = [  9.495700d+02,  8.364100d-03,  1.398800d-04,  1.498700d+00, 0.d0, 0.d0]
    self%cls(024,:) = [  9.639900d+02,  9.726500d-03,  1.508800d-04,  1.559300d+00, 0.d0, 0.d0]
    self%cls(025,:) = [  9.787400d+02,  1.126100d-02,  1.622300d-04,  1.617700d+00, 0.d0, 0.d0]
    self%cls(026,:) = [  9.938100d+02,  1.297000d-02,  1.739200d-04,  1.673100d+00, 0.d0, 0.d0]
    self%cls(027,:) = [  1.009200d+03,  1.485600d-02,  1.859700d-04,  1.724700d+00, 0.d0, 0.d0]
    self%cls(028,:) = [  1.024900d+03,  1.692100d-02,  1.983900d-04,  1.771800d+00, 0.d0, 0.d0]
    self%cls(029,:) = [  1.040800d+03,  1.916800d-02,  2.111700d-04,  1.813400d+00, 0.d0, 0.d0]
    self%cls(030,:) = [  1.057100d+03,  2.159900d-02,  2.243500d-04,  1.848800d+00, 0.d0, 0.d0]
    self%cls(031,:) = [  1.073600d+03,  2.421700d-02,  2.379300d-04,  1.877300d+00, 0.d0, 0.d0]
    self%cls(032,:) = [  1.090300d+03,  2.702800d-02,  2.519300d-04,  1.898200d+00, 0.d0, 0.d0]
    self%cls(033,:) = [  1.107300d+03,  3.003900d-02,  2.663800d-04,  1.911000d+00, 0.d0, 0.d0]
    self%cls(034,:) = [  1.124600d+03,  3.325500d-02,  2.812900d-04,  1.915300d+00, 0.d0, 0.d0]
    self%cls(035,:) = [  1.142000d+03,  3.668500d-02,  2.967000d-04,  1.910400d+00, 0.d0, 0.d0]
    self%cls(036,:) = [  1.159700d+03,  4.033400d-02,  3.126300d-04,  1.895900d+00, 0.d0, 0.d0]
    self%cls(037,:) = [  1.177700d+03,  4.421000d-02,  3.291000d-04,  1.871200d+00, 0.d0, 0.d0]
    self%cls(038,:) = [  1.195800d+03,  4.831800d-02,  3.461500d-04,  1.835900d+00, 0.d0, 0.d0]
    self%cls(039,:) = [  1.214100d+03,  5.266500d-02,  3.638100d-04,  1.789200d+00, 0.d0, 0.d0]
    self%cls(040,:) = [  1.232600d+03,  5.725900d-02,  3.821000d-04,  1.730900d+00, 0.d0, 0.d0]
    self%cls(041,:) = [  1.251300d+03,  6.210500d-02,  4.010700d-04,  1.660300d+00, 0.d0, 0.d0]
    self%cls(042,:) = [  1.270200d+03,  6.720800d-02,  4.207300d-04,  1.577400d+00, 0.d0, 0.d0]
    self%cls(043,:) = [  1.289200d+03,  7.257200d-02,  4.411200d-04,  1.482200d+00, 0.d0, 0.d0]
    self%cls(044,:) = [  1.308500d+03,  7.820200d-02,  4.622600d-04,  1.374800d+00, 0.d0, 0.d0]
    self%cls(045,:) = [  1.327900d+03,  8.410200d-02,  4.842000d-04,  1.255000d+00, 0.d0, 0.d0]
    self%cls(046,:) = [  1.347500d+03,  9.027700d-02,  5.069500d-04,  1.122900d+00, 0.d0, 0.d0]
    self%cls(047,:) = [  1.367400d+03,  9.673100d-02,  5.305400d-04,  9.785200d-01, 0.d0, 0.d0]
    self%cls(048,:) = [  1.387400d+03,  1.034700d-01,  5.549900d-04,  8.218100d-01, 0.d0, 0.d0]
    self%cls(049,:) = [  1.407600d+03,  1.104900d-01,  5.803200d-04,  6.527900d-01, 0.d0, 0.d0]
    self%cls(050,:) = [  1.428100d+03,  1.178100d-01,  6.065500d-04,  4.714400d-01, 0.d0, 0.d0]
    self%cls(051,:) = [  1.448700d+03,  1.254200d-01,  6.337000d-04,  2.777400d-01, 0.d0, 0.d0]
    self%cls(052,:) = [  1.469600d+03,  1.333300d-01,  6.617600d-04,  7.160200d-02, 0.d0, 0.d0]
    self%cls(053,:) = [  1.490600d+03,  1.415300d-01,  6.907700d-04, -1.470900d-01, 0.d0, 0.d0]
    self%cls(054,:) = [  1.511900d+03,  1.500100d-01,  7.207100d-04, -3.784600d-01, 0.d0, 0.d0]
    self%cls(055,:) = [  1.533400d+03,  1.587900d-01,  7.515800d-04, -6.226100d-01, 0.d0, 0.d0]
    self%cls(056,:) = [  1.555000d+03,  1.678500d-01,  7.833900d-04, -8.796600d-01, 0.d0, 0.d0]
    self%cls(057,:) = [  1.576900d+03,  1.771900d-01,  8.161200d-04, -1.149700d+00, 0.d0, 0.d0]
    self%cls(058,:) = [  1.599000d+03,  1.868100d-01,  8.497700d-04, -1.432900d+00, 0.d0, 0.d0]
    self%cls(059,:) = [  1.621300d+03,  1.967200d-01,  8.843100d-04, -1.729300d+00, 0.d0, 0.d0]
    self%cls(060,:) = [  1.643800d+03,  2.068900d-01,  9.197200d-04, -2.039100d+00, 0.d0, 0.d0]
    self%cls(061,:) = [  1.666500d+03,  2.173400d-01,  9.560000d-04, -2.362400d+00, 0.d0, 0.d0]
    self%cls(062,:) = [  1.689300d+03,  2.280600d-01,  9.930900d-04, -2.699000d+00, 0.d0, 0.d0]
    self%cls(063,:) = [  1.712400d+03,  2.390500d-01,  1.031000d-03, -3.049000d+00, 0.d0, 0.d0]
    self%cls(064,:) = [  1.735800d+03,  2.503100d-01,  1.069600d-03, -3.412200d+00, 0.d0, 0.d0]
    self%cls(065,:) = [  1.759300d+03,  2.618400d-01,  1.109000d-03, -3.788600d+00, 0.d0, 0.d0]
    self%cls(066,:) = [  1.783000d+03,  2.736300d-01,  1.149100d-03, -4.178200d+00, 0.d0, 0.d0]
    self%cls(067,:) = [  1.807000d+03,  2.856800d-01,  1.189800d-03, -4.580800d+00, 0.d0, 0.d0]
    self%cls(068,:) = [  1.831200d+03,  2.980000d-01,  1.231000d-03, -4.996300d+00, 0.d0, 0.d0]
    self%cls(069,:) = [  1.855600d+03,  3.105800d-01,  1.272900d-03, -5.424800d+00, 0.d0, 0.d0]
    self%cls(070,:) = [  1.880200d+03,  3.234100d-01,  1.315200d-03, -5.866200d+00, 0.d0, 0.d0]
    self%cls(071,:) = [  1.905100d+03,  3.365000d-01,  1.358000d-03, -6.320200d+00, 0.d0, 0.d0]
    self%cls(072,:) = [  1.930200d+03,  3.498300d-01,  1.401200d-03, -6.786800d+00, 0.d0, 0.d0]
    self%cls(073,:) = [  1.955600d+03,  3.634000d-01,  1.444800d-03, -7.265400d+00, 0.d0, 0.d0]
    self%cls(074,:) = [  1.981200d+03,  3.772000d-01,  1.488600d-03, -7.755800d+00, 0.d0, 0.d0]
    self%cls(075,:) = [  2.007100d+03,  3.912100d-01,  1.532700d-03, -8.257700d+00, 0.d0, 0.d0]
    self%cls(076,:) = [  2.033100d+03,  4.054300d-01,  1.577100d-03, -8.770800d+00, 0.d0, 0.d0]
    self%cls(077,:) = [  2.059500d+03,  4.198500d-01,  1.621600d-03, -9.294600d+00, 0.d0, 0.d0]
    self%cls(078,:) = [  2.086000d+03,  4.344500d-01,  1.666300d-03, -9.828900d+00, 0.d0, 0.d0]
    self%cls(079,:) = [  2.112800d+03,  4.492300d-01,  1.711000d-03, -1.037300d+01, 0.d0, 0.d0]
    self%cls(080,:) = [  2.139900d+03,  4.641700d-01,  1.755800d-03, -1.092800d+01, 0.d0, 0.d0]
    self%cls(081,:) = [  2.167200d+03,  4.792600d-01,  1.800700d-03, -1.149100d+01, 0.d0, 0.d0]
    self%cls(082,:) = [  2.194700d+03,  4.944900d-01,  1.845500d-03, -1.206400d+01, 0.d0, 0.d0]
    self%cls(083,:) = [  2.222500d+03,  5.098400d-01,  1.890400d-03, -1.264500d+01, 0.d0, 0.d0]
    self%cls(084,:) = [  2.250500d+03,  5.252700d-01,  1.935200d-03, -1.323500d+01, 0.d0, 0.d0]
    self%cls(085,:) = [  2.278800d+03,  5.407800d-01,  1.980000d-03, -1.383200d+01, 0.d0, 0.d0]
    self%cls(086,:) = [  2.307300d+03,  5.563400d-01,  2.024800d-03, -1.443600d+01, 0.d0, 0.d0]
    self%cls(087,:) = [  2.336100d+03,  5.719300d-01,  2.069500d-03, -1.504700d+01, 0.d0, 0.d0]
    self%cls(088,:) = [  2.365100d+03,  5.875300d-01,  2.114100d-03, -1.566500d+01, 0.d0, 0.d0]
    self%cls(089,:) = [  2.394300d+03,  6.031100d-01,  2.158800d-03, -1.628800d+01, 0.d0, 0.d0]
    self%cls(090,:) = [  2.423800d+03,  6.186600d-01,  2.203400d-03, -1.691600d+01, 0.d0, 0.d0]
    self%cls(091,:) = [  2.453500d+03,  6.341600d-01,  2.248100d-03, -1.754900d+01, 0.d0, 0.d0]
    self%cls(092,:) = [  2.483400d+03,  6.495900d-01,  2.292800d-03, -1.818700d+01, 0.d0, 0.d0]
    self%cls(093,:) = [  2.513600d+03,  6.649500d-01,  2.337500d-03, -1.882900d+01, 0.d0, 0.d0]
    self%cls(094,:) = [  2.544000d+03,  6.802200d-01,  2.382400d-03, -1.947400d+01, 0.d0, 0.d0]
    self%cls(095,:) = [  2.574600d+03,  6.953900d-01,  2.427500d-03, -2.012300d+01, 0.d0, 0.d0]
    self%cls(096,:) = [  2.605500d+03,  7.104600d-01,  2.472800d-03, -2.077400d+01, 0.d0, 0.d0]
    self%cls(097,:) = [  2.636600d+03,  7.254100d-01,  2.518300d-03, -2.142700d+01, 0.d0, 0.d0]
    self%cls(098,:) = [  2.667800d+03,  7.402300d-01,  2.564100d-03, -2.208100d+01, 0.d0, 0.d0]
    self%cls(099,:) = [  2.699300d+03,  7.549100d-01,  2.610400d-03, -2.273800d+01, 0.d0, 0.d0]
    self%cls(100,:) = [  2.731000d+03,  7.694300d-01,  2.657000d-03, -2.339400d+01, 0.d0, 0.d0]
    self%cls(101,:) = [  2.762900d+03,  7.838000d-01,  2.704200d-03, -2.405200d+01, 0.d0, 0.d0]
    self%cls(102,:) = [  2.795000d+03,  7.980000d-01,  2.751900d-03, -2.470900d+01, 0.d0, 0.d0]
    self%cls(103,:) = [  2.827200d+03,  8.120100d-01,  2.800300d-03, -2.536500d+01, 0.d0, 0.d0]
    self%cls(104,:) = [  2.859700d+03,  8.258300d-01,  2.849400d-03, -2.602100d+01, 0.d0, 0.d0]
    self%cls(105,:) = [  2.892300d+03,  8.394400d-01,  2.899300d-03, -2.667500d+01, 0.d0, 0.d0]
    self%cls(106,:) = [  2.925100d+03,  8.528400d-01,  2.950000d-03, -2.732700d+01, 0.d0, 0.d0]
    self%cls(107,:) = [  2.958100d+03,  8.660200d-01,  3.001600d-03, -2.797600d+01, 0.d0, 0.d0]
    self%cls(108,:) = [  2.991300d+03,  8.789500d-01,  3.054200d-03, -2.862300d+01, 0.d0, 0.d0]
    self%cls(109,:) = [  3.024600d+03,  8.916500d-01,  3.107800d-03, -2.926600d+01, 0.d0, 0.d0]
    self%cls(110,:) = [  3.058100d+03,  9.040800d-01,  3.162600d-03, -2.990600d+01, 0.d0, 0.d0]
    self%cls(111,:) = [  3.091700d+03,  9.162400d-01,  3.218500d-03, -3.054200d+01, 0.d0, 0.d0]
    self%cls(112,:) = [  3.125500d+03,  9.281200d-01,  3.275600d-03, -3.117200d+01, 0.d0, 0.d0]
    self%cls(113,:) = [  3.159400d+03,  9.397100d-01,  3.333900d-03, -3.179700d+01, 0.d0, 0.d0]
    self%cls(114,:) = [  3.193500d+03,  9.510000d-01,  3.393600d-03, -3.241700d+01, 0.d0, 0.d0]
    self%cls(115,:) = [  3.227600d+03,  9.619800d-01,  3.454600d-03, -3.302900d+01, 0.d0, 0.d0]
    self%cls(116,:) = [  3.261900d+03,  9.726300d-01,  3.516900d-03, -3.363500d+01, 0.d0, 0.d0]
    self%cls(117,:) = [  3.296300d+03,  9.829500d-01,  3.580600d-03, -3.423200d+01, 0.d0, 0.d0]
    self%cls(118,:) = [  3.330800d+03,  9.929200d-01,  3.645800d-03, -3.482200d+01, 0.d0, 0.d0]
    self%cls(119,:) = [  3.365400d+03,  1.002500d+00,  3.712300d-03, -3.540200d+01, 0.d0, 0.d0]
    self%cls(120,:) = [  3.400100d+03,  1.011800d+00,  3.780300d-03, -3.597300d+01, 0.d0, 0.d0]
    self%cls(121,:) = [  3.434800d+03,  1.020600d+00,  3.849700d-03, -3.653300d+01, 0.d0, 0.d0]
    self%cls(122,:) = [  3.469600d+03,  1.029100d+00,  3.920500d-03, -3.708300d+01, 0.d0, 0.d0]
    self%cls(123,:) = [  3.504500d+03,  1.037200d+00,  3.992700d-03, -3.762100d+01, 0.d0, 0.d0]
    self%cls(124,:) = [  3.539400d+03,  1.044800d+00,  4.066300d-03, -3.814800d+01, 0.d0, 0.d0]
    self%cls(125,:) = [  3.574300d+03,  1.052100d+00,  4.141300d-03, -3.866100d+01, 0.d0, 0.d0]
    self%cls(126,:) = [  3.609300d+03,  1.058800d+00,  4.217500d-03, -3.916200d+01, 0.d0, 0.d0]
    self%cls(127,:) = [  3.644300d+03,  1.065200d+00,  4.295000d-03, -3.964900d+01, 0.d0, 0.d0]
    self%cls(128,:) = [  3.679300d+03,  1.071100d+00,  4.373800d-03, -4.012200d+01, 0.d0, 0.d0]
    self%cls(129,:) = [  3.714300d+03,  1.076500d+00,  4.453700d-03, -4.057900d+01, 0.d0, 0.d0]
    self%cls(130,:) = [  3.749300d+03,  1.081400d+00,  4.534800d-03, -4.102100d+01, 0.d0, 0.d0]
  end subroutine init

end module fiducial_cls_mod_smw

module bflike_utils_smw

  use healpix_types
  implicit none

  private
  integer(i4b) ,parameter ,public :: smwTT=1 ,smwEE=2 ,smwBB=3 ,smwTE=4 ,smwTB=5 ,smwEB=6

  type allocatable_vector
     real(dp) ,allocatable :: row(:)
  end type allocatable_vector
  
  type triangular_matrix

     integer :: size 
     logical :: allocated 
     type(allocatable_vector) ,allocatable :: column(:)

   contains
     procedure :: alloc => alloc_tm
     procedure :: clean => clean_tm

  end type triangular_matrix

  type allocatable_matrix

     integer  :: nrow
     integer  :: ncol
     logical  :: allocated
     real(dp) ,allocatable :: m(:,:)

   contains

     procedure :: alloc => alloc_am
     procedure :: clean => clean_am

  end type allocatable_matrix


  public allocatable_matrix ,triangular_matrix

contains

!allocation and deallocation procedures for derived types
  subroutine alloc_tm(self,n)
    class(triangular_matrix) ,intent(inout) :: self
    integer ,intent(in)                     :: n
    
    integer                                 :: i
    
    if(self%allocated) call self%clean()
    
    self%size = n
    allocate(self%column(n))
    
    do i = 1,n
       allocate(self%column(i)%row(i:n))
    end do
    
    self%allocated = .true.
    
  end subroutine alloc_tm

  subroutine clean_tm(self)
    class(triangular_matrix) ,intent(inout) :: self
    integer                                 :: i
    if(.not. self%allocated) return
    
    do i = 1,self%size
       deallocate(self%column(i)%row)
    end do
    deallocate(self%column)
    
    self%size = -1
    self%allocated = .false.
    
  end subroutine clean_tm

  subroutine alloc_am(self,nrow,ncol)
    class(allocatable_matrix) ,intent(inout) :: self
    integer ,intent(in)                     :: nrow,ncol
      
    if(self%allocated) call self%clean()
    
    self%nrow = nrow
    self%nrow = ncol
    allocate(self%m(nrow,ncol))
    
    self%allocated = .true.
    
  end subroutine alloc_am

  subroutine clean_am(self)
    class(allocatable_matrix) ,intent(inout) :: self
    if(.not. self%allocated) return

    deallocate(self%m)
    
    self%nrow = -1
    self%ncol = -1
    self%allocated = .false.
    
  end subroutine clean_am

end module bflike_utils_smw

module bflike_smw

  use healpix_types
  use bflike_utils_smw
  use fitstools_smw , only : read_dbintab ,printerror 
  implicit none

  private

  public init_pix_like_smw ,clean_pix_like_smw ,get_pix_loglike_smw ,external_compute_plms

  integer,save :: ntemp ,nq ,nu ,ntot ,nqu ,lmax ,lswitch ,ndata

  real(dp) ,allocatable ,dimension(:,:) ,save :: clnorm ,Tvec ,Qvec ,Uvec ,dt ,auxdt ,cls
  real(dp) ,allocatable ,dimension(:) ,save   :: reflike ,pl,plm,f1,f2
  real(dp) ,allocatable ,target               :: S(:,:)

  type(triangular_matrix) ,save               :: cos1 ,cos2 ,sin1 ,sin2 ,ncov 
  type(allocatable_matrix) ,allocatable ,save :: fevec(:) ,feval(:) ,cblocks(:)

  abstract interface

     subroutine get_cov_interface(clsin,linf,lsup,cov,project_mondip,symmetrize)

       import dp
       implicit none
       real(dp) ,intent(in) :: clsin(2:,:)
       integer ,intent(in)  :: linf ,lsup
       real(dp) ,intent(out):: cov(:,:)
       logical ,optional ,intent(in) :: project_mondip ,symmetrize
       
     end subroutine get_cov_interface

  end  interface
    
contains

  subroutine external_compute_plms(passntemp,passnq,passnu,theta,phi,beam,passlswitch,basisfile)
    implicit none
    integer ,intent(in)  :: passntemp ,passnq ,passnu ,passlswitch
    real(dp) ,intent(in) :: theta(:) ,phi(:) ,beam(0:,:) 
    character(len=*) ,intent(in) :: basisfile
    
    real(dp)            :: mod
    integer             :: i ,j ,jdum ,l

    ntemp = passntemp
    nq    = passnq
    nu    = passnu
    lswitch = passlswitch
    lmax    = lswitch

    ntot = ntemp +nq +nu
    nqu  = nq +nu

    allocate(Tvec(3,ntemp),Qvec(3,nq),Uvec(3,nu))
    jdum = -1
    do i=1,ntemp
       Tvec(1,i) = sin(theta(i))*cos(phi(i))
       Tvec(2,i) = sin(theta(i))*sin(phi(i))
       Tvec(3,i) = cos(theta(i))
       mod = sqrt(sum(Tvec(:,i)**2))
       Tvec(:,i) = Tvec(:,i)/mod
    end do
    do i=1,nq
       Qvec(1,i) = sin(theta(i+ntemp))*cos(phi(i+ntemp))
       Qvec(2,i) = sin(theta(i+ntemp))*sin(phi(i+ntemp))
       Qvec(3,i) = cos(theta(i+ntemp))
       mod = sqrt(sum(Qvec(:,i)**2))
       Qvec(:,i) = Qvec(:,i)/mod
    end do
    do i=1,nu
       Uvec(1,i) = sin(theta(i+ntemp+nq))*cos(phi(i+ntemp+nq))
       Uvec(2,i) = sin(theta(i+ntemp+nq))*sin(phi(i+ntemp+nq))
       Uvec(3,i) = cos(theta(i+ntemp+nq))
       mod = sqrt(sum(Uvec(:,i)**2))
       Uvec(:,i) = Uvec(:,i)/mod
    end do

    allocate(clnorm(2:lmax,6))
    call set_clnorm(beam(2:,:),clnorm(2:,:))

    call cos1%alloc(ntot)
    call cos2%alloc(ntot)
    call sin1%alloc(ntot)
    call sin2%alloc(ntot)

    allocate(cls(2:lmax,6))
    allocate(pl(1:lmax),plm(1:lmax))
    allocate(f1(2:lmax),f2(2:lmax))

    call precompute_rotation_angle()

    allocate(fevec(2:lswitch),feval(2:lswitch),cblocks(2:lswitch))
    do l=2,lswitch
       call fevec(l)%alloc(nrow=ntot,ncol=3*(2*l+1))
       fevec(l)%m = 0.d0
       call feval(l)%alloc(nrow=3*(2*l+1),ncol=1)
       feval(l)%m = 0.d0
       call cblocks(l)%alloc(nrow=3*(2*l+1),ncol=3*(2*l+1))
       cblocks(l)%m = 0.d0
    end do

    if(.not.read_plms_from_file(trim(basisfile),feval,fevec,cblocks)) then
       write(*,*) 'basis file not found/wrong number of modes'
       write(*,*) 'computing plms from scratch'
       write(*,*) 'computed plms will be saved in: '//trim(basisfile)
       call compute_plms(feval,fevec,cblocks,trim(basisfile))
       write(*,*) 'done'
    end if

    deallocate(Tvec,Qvec,Uvec)
    call cos1%clean()
    call cos2%clean()
    call sin1%clean()
    call sin2%clean()

    do l=2,lswitch
       call fevec(l)%clean()
       call feval(l)%clean()
       call cblocks(l)%clean()
    end do
    
    deallocate(fevec,feval,cblocks)

    deallocate(cls,pl,plm,f1,f2)

  end subroutine external_compute_plms

  subroutine init_pix_like_smw(clik_bflike_dir)

    implicit none
    character(len=*) ,intent(in) :: clik_bflike_dir
    
    integer(i4b)          :: i ,j ,l ,myunit ,err  &
         ,iq ,iu ,jq ,ju ,il ,nside ,nwrite ,jdum ,blocksize ,readwrite ,hdutype
    real(dp)              :: theta ,phi ,ell ,mod ,a1 ,a2 ,nullval
    real(dp) ,allocatable :: evec(:,:) ,clstmp(:,:) ,NCVM(:,:) ,bl(:,:) 
    character(len=256)    :: datafile ,basisfile ,clfiducial
    character(len=80)     :: comment
    logical               :: project_mondip ,anynull
    integer(i4b)          :: info ,neigen 
    
    namelist /inputs/ datafile ,project_mondip ,lmax ,lswitch ,basisfile ,clfiducial 

    myunit = 42

    open(myunit,file=trim(clik_bflike_dir)//'/params_bflike.ini',status='old')
    read(myunit,inputs)
    close(myunit)

    err = 0
    readwrite = 0

    call ftopen(myunit,trim(clik_bflike_dir)//datafile,readwrite,blocksize,err)
    if(err .gt. 0) then 
       write(*,*) 'OPEN'
       call printerror(info)
       stop
    end if

    call ftmrhd(myunit, +1, hdutype, err)
    if(err >0) then 
       write(*,*) 'HDR'
       call printerror(err)
       stop
    end if

    call ftgkyj(myunit,'NTEMP',ntemp ,comment ,err)
    if(err >0) then
       write(*,*) 'NTEMP'
       call printerror(err)
       stop
    end if

    call ftgkyj(myunit,'NQ',nq ,comment ,err)
    if(err >0) then
       write(*,*) 'NQ'
       call printerror(err)
       stop
    end if
    call ftgkyj(myunit,'NU',nu ,comment ,err)
    if(err >0) then
       write(*,*) 'NU'
       call printerror(err)
       stop
    end if

    ntot = ntemp +nq +nu
    nqu  = nq +nu
    
    call ftgkyj(myunit,'NSIDE',nside ,comment ,err)
    if(err >0) then 
       write(*,*) 'NSIDE'
       call printerror(err)
       stop
    end if
    if(lmax .gt. 4*nside) then
       write(*,*) 'WARNING: you have nside =',nside,' and lmax =',lmax
       write(*,*) 'WARNING: Setting lmax =4*nside'
       lmax = 4*nside
    end if

    call ftgkyj(myunit,'NUMDATA',ndata ,comment ,err)
    if(err >0) then
       write(*,*) 'NUMDATA'
       call printerror(err)
       stop
    end if
    if(ndata .gt. 1) then
       write(*,'(a)') 'WARNING: This version allows only 1 map '
       write(*,'(a)') 'WARNING: Setting numdata = 1'
       ndata = 1
    end if

    call ftgkyj(myunit,'NWRITE',nwrite ,comment ,err)
    if(err >0) then
       write(*,*) 'NWRITE'
       call printerror(err)
       stop
    end if

    call ftclos(myunit, err)
    if (err > 0) call printerror(err)
    
    write(*,*) 'BFLike Ntemp  =',ntemp
    write(*,*) 'BFLike Nq     =',nq
    write(*,*) 'BFLike Nu     =',nu
    write(*,*) 'BFLike Nside  =',nside
    write(*,*) 'BFLike Nwrite =',nwrite

    allocate(evec(0:nwrite-1,1))
    call read_dbintab(trim(clik_bflike_dir)//datafile, evec, nwrite, 1, nullval, anynull)
    if(anynull) then
       write(*,*) 'Missing pixels have been given the value ',nullval
    end if

!geometry
    allocate(Tvec(3,ntemp),Qvec(3,nq),Uvec(3,nu))
    jdum = -1
    do i=1,ntemp
       jdum = jdum+1;       theta = evec(jdum,1)
       jdum = jdum+1;       phi = evec(jdum,1)
       Tvec(1,i) = sin(theta)*cos(phi)
       Tvec(2,i) = sin(theta)*sin(phi)
       Tvec(3,i) = cos(theta)
       mod = sqrt(sum(Tvec(:,i)**2))
       Tvec(:,i) = Tvec(:,i)/mod
    end do
    do i=1,nq
       jdum = jdum+1;       theta = evec(jdum,1)
       jdum = jdum+1;       phi = evec(jdum,1)
       Qvec(1,i) = sin(theta)*cos(phi)
       Qvec(2,i) = sin(theta)*sin(phi)
       Qvec(3,i) = cos(theta)
       mod = sqrt(sum(Qvec(:,i)**2))
       Qvec(:,i) = Qvec(:,i)/mod
    end do
    do i=1,nu
       jdum = jdum+1;       theta = evec(jdum,1)
       jdum = jdum+1;       phi = evec(jdum,1)
       Uvec(1,i) = sin(theta)*cos(phi)
       Uvec(2,i) = sin(theta)*sin(phi)
       Uvec(3,i) = cos(theta)
       mod = sqrt(sum(Uvec(:,i)**2))
       Uvec(:,i) = Uvec(:,i)/mod
    end do
       
!beam
    iu = jdum
    allocate(bl(0:4*nside,6),clnorm(2:lmax,6))
    do i=1,6
       il = iu+1
       iu = il +4*nside
       bl(0:4*nside,i) = evec(il:iu,1)
    end do

    call set_clnorm(bl(2:,:),clnorm(2:,:))
    deallocate(bl)

!ncvm
    allocate(ncvm(ntot,ntot))
    do i=1,ntot !ncm
       il = iu +1
       iu = il +ntot-1
       ncvm(:,i) = evec(il:iu,1)
    end do

!data
    allocate(dt(ntot,ndata))
    do j=1,ndata
       il = iu +1
       iu = il +ntot-1
       dt(:,1) = evec(il:iu,1)
    end do
    deallocate(evec)

    allocate(reflike(ndata))

    allocate(cls(2:lmax,6))
    allocate(pl(1:lmax),plm(1:lmax))
    allocate(f1(2:lmax),f2(2:lmax))
    allocate(auxdt(ntot,ndata))

!precompute the rotation angles
    call cos1%alloc(ntot)
    call cos2%alloc(ntot)
    call sin1%alloc(ntot)
    call sin2%alloc(ntot)

    call precompute_rotation_angle()

!update the NCVM
    allocate(clstmp(2:lmax,6))
    clstmp = 0._dp
    
    if(lswitch .lt. lmax) then
       if(.not. read_cls_from_file(clfiducial,clstmp)) then
          write(0,*) 'WARNING: '//trim(clfiducial)//' not found or not enough columns'
          write(0,*) 'using default values'
       end if
       call update_ncvm(clstmp,lswitch+1,lmax,NCVM,project_mondip=project_mondip)
    else
       call update_ncvm(clstmp,2,2,NCVM,project_mondip=project_mondip)
    end if
    deallocate(clstmp)

!
    if(3*(lswitch+1)**2 .gt. ntot) then
       write(*,*) 'requested number of modes is higher than the total number of pixels'
       write(*,*) 'try lowering lswitch. usually 3*(lswitch+1)^2 ~ 0.5-0.8 n_pix is a good compromise'
       stop
    end if

    allocate(fevec(2:lswitch),feval(2:lswitch),cblocks(2:lswitch))
    do l=2,lswitch
       call fevec(l)%alloc(nrow=ntot,ncol=3*(2*l+1))
       fevec(l)%m = 0.d0
       call feval(l)%alloc(nrow=3*(2*l+1),ncol=1)
       feval(l)%m = 0.d0
       call cblocks(l)%alloc(nrow=3*(2*l+1),ncol=3*(2*l+1))
       cblocks(l)%m = 0.d0
    end do

    !look for plms eigenmodes
    if(.not.read_plms_from_file(basisfile,feval,fevec,cblocks)) then
       write(*,*) 'basis file not found/wrong number of modes'
       write(*,*) 'computing plms from scratch'
       write(*,*) 'computed plms will be saved in: '//trim(basisfile)
       call compute_plms(feval,fevec,cblocks,trim(basisfile))
       write(*,*) 'done'
    end if

    !contract NCVM and data on the plms 
    neigen = 0
    do l=2,lswitch
       neigen = neigen +2*l+1
    end do
    neigen = neigen*3

    allocate(evec(ntot,neigen))
    iu = 0
    do l=2,lswitch
       il = iu +1
       iu = il +3*(2*l +1)-1
       evec(:,il:iu) = fevec(l)%m
    end do
    !!!allocate(S,source=evec)
    allocate(S(lbound(evec,1):ubound(evec,1), lbound(evec,2):ubound(evec,2)), source=evec)

    auxdt = dt
    call dposv('L',ntot,ndata,NCVM,ntot,auxdt,ntot,info)
    if(info .ne. 0) then 
       write(*,*) '(updated) NCVM is not invertible'
       write(*,*) 'info = ',info
       stop
    end if

    do i = 1,ndata
       reflike(i) = sum(dt(:,i)*auxdt(:,i))
    end do

    deallocate(dt)  ;allocate(dt(neigen,ndata))
    call dgemm('T','N', neigen, ndata, ntot, 1.d0, evec, &
         ntot, auxdt, ntot, 0.d0, dt, neigen)
    !deallocate(auxdt) ;allocate(auxdt,source = dt)
    deallocate(auxdt)
    allocate(auxdt(lbound(dt,1):ubound(dt,1), lbound(dt,2):ubound(dt,2)), source=dt)

    call dpotrs('L',ntot,neigen,NCVM,ntot,evec,ntot,info)
    write(*,*) 'info = ',info

    deallocate(NCVM) ;allocate(NCVM(neigen,neigen))
    call dgemm('T','N', neigen, neigen, ntot, 1.d0, S, &
         ntot, evec, ntot, 0.d0, NCVM , neigen)

    call ncov%alloc(neigen)
    do i = 1,neigen
       ncov%column(i)%row(i:neigen) = ncvm(i:neigen,i)
    end do
    deallocate(NCVM)

    deallocate(evec)
    deallocate(S) ;allocate(S(neigen,neigen))
    
    !clean up all the stuff we don't need anymore
    do l=2,lswitch
       call fevec(l)%clean()
    end do
    deallocate(fevec)

    deallocate(Tvec,Uvec,Qvec)
    deallocate(pl,plm,f1,f2)

    call cos1%clean()
    call cos2%clean()
    call sin1%clean()
    call sin2%clean()
    deallocate(clnorm)

  end subroutine init_pix_like_smw

  subroutine precompute_rotation_angle()
    implicit none

    integer(i4b) :: i ,j ,iq ,jq ,iu ,ju
    real(dp)     :: a1 ,a2

! QT
    do i=1,ntemp
       do j=1,nq
          call get_rotation_angle(QVec(:,j),TVec(:,i),a1,a2)
          jq = j +ntemp
          cos1%column(i)%row(jq) = cos(a1)
          sin1%column(i)%row(jq) = sin(a1)

       enddo
!UT
       do j=1,nu
          call get_rotation_angle(UVec(:,j),TVec(:,i),a1,a2)
          ju = j +ntemp +nq
          cos1%column(i)%row(ju) = cos(a1)
          sin1%column(i)%row(ju) = sin(a1)
       end do
    end do

    do i = 1,nq
!QQ
       iq = i+ntemp
       do j = i,nq
          call get_rotation_angle(QVec(:,j),QVec(:,i),a1,a2)
          jq = j+ntemp

          cos1%column(iq)%row(jq) = cos(a1)
          sin1%column(iq)%row(jq) = sin(a1)
          cos2%column(iq)%row(jq) = cos(a2)
          sin2%column(iq)%row(jq) = sin(a2)

       end do

       do j=1,nu
          call get_rotation_angle(UVec(:,j),QVec(:,i),a1,a2)

          ju = j+ntemp+nq

          cos1%column(iq)%row(ju) = cos(a1)
          sin1%column(iq)%row(ju) = sin(a1)
          cos2%column(iq)%row(ju) = cos(a2)
          sin2%column(iq)%row(ju) = sin(a2)

       end do
    end do

    do i =1,nu
!UU
       iu = i+ntemp+nq
       do j=i,nu

          call get_rotation_angle(UVec(:,j),UVec(:,i),a1,a2)

          ju = j+ntemp+nq

          cos1%column(iu)%row(ju) = cos(a1)
          sin1%column(iu)%row(ju) = sin(a1)
          cos2%column(iu)%row(ju) = cos(a2)
          sin2%column(iu)%row(ju) = sin(a2)

       end do
    end do

  end subroutine precompute_rotation_angle

  subroutine set_clnorm(bl,clnorm)
    implicit none
    real(dp) ,intent(in)  :: bl(2:,:)
    real(dp) ,intent(out) :: clnorm(2:,:)
    
    integer               :: l
    real(dp)              :: fct ,fct2 ,ell ,chngconv

    fct2 = 1._dp/24._dp
    fct = sqrt(fct2)
    
    clnorm(2,smwTT) = bl(2,smwTT)**2/2.4_dp
    clnorm(2,smwEE) = bl(2,smwEE)**2*fct2/2.4_dp
    clnorm(2,smwBB) = bl(2,smwBB)**2*fct2/2.4_dp
    clnorm(2,smwTE) = bl(2,smwTE)**2*fct/2.4_dp
    clnorm(2,smwTB) = bl(2,smwTB)**2*fct/2.4_dp
    clnorm(2,smwEB) = bl(2,smwEB)**2*fct2/2.4_dp

    do l=3,lmax
       ell = real(l,kind=dp)

       fct2 = 1.d0/((ell+2.d0)*(ell+1.d0)*ell*(ell-1.d0))
       fct = sqrt(fct2)
       chngconv = (2.d0*ell +1.d0)/2.d0/(ell*(ell+1.d0))

       clnorm(l,smwTT) = chngconv*bl(l,smwTT)**2
       clnorm(l,smwEE) = chngconv*bl(l,smwEE)**2*fct2
       clnorm(l,smwBB) = chngconv*bl(l,smwBB)**2*fct2
       clnorm(l,smwTE) = chngconv*bl(l,smwTE)**2*fct
       clnorm(l,smwTB) = chngconv*bl(l,smwTB)**2*fct
       clnorm(l,smwEB) = chngconv*bl(l,smwEB)**2*fct2

    end do

  end subroutine set_clnorm

  subroutine compute_plms(eigenval,eigenvec,blocks,file)

    implicit none
    type(allocatable_matrix) ,intent(inout) :: eigenval(2:) ,eigenvec(2:) ,blocks(2:)
    character(len=*) ,intent(in) ,optional :: file

    integer(i4b)              :: il ,iu ,l ,lwork ,liwork ,info ,neigenl 
    real(dp)     ,allocatable :: work(:) ,evec(:,:) ,eval(:) ,clstmp(:,:) &
         ,Z(:,:)
    integer(i4b) ,allocatable :: isupp(:) ,iwork(:)
 
    double precision dlamch
    external dlamch

    allocate(Z(ntot,ntot))

    !TT plms and eigenvalues
    allocate(evec(ntemp,ntemp),eval(ntemp),isupp(2*ntemp),work(1),iwork(1))
    call dsyevr('V', 'A', 'L', ntemp, Z(1:ntemp,1:ntemp), ntemp, 0.d0, 0.d0, 0, 0, dlamch('S'), neigenl, eval, evec, ntemp, isupp, work, -1, iwork, -1, info)
    lwork = work(1) ;liwork = iwork(1)
    deallocate(work,iwork) ;allocate(work(lwork),iwork(liwork))

    allocate(clstmp(2:lmax,6))
    Z = 0.d0
    do l =2,lswitch
       write(*,*) 'computing basis for TT l=',l
       il = ntemp +1 -(2*l+1)
       iu = ntemp
       
       clstmp = 0._dp
       clstmp(l,smwTT) = 1._dp
       call get_tt_cov(clstmp,l,l,Z(1:ntemp,1:ntemp),project_mondip=.false.)
       call dsyevr('V', 'I', 'L', ntemp, Z(1:ntemp,1:ntemp), ntemp, 0.d0, 0.d0, il, iu, 2*dlamch('S'), neigenl, eval, evec, ntemp, isupp, work, lwork, iwork, liwork, info)
       eigenvec(l)%m(1:ntemp,1:neigenl) = evec(1:ntemp,1:neigenl)
       eigenval(l)%m(1:neigenl,1)       = eval(1:neigenl)

    end do
    deallocate(evec,eval,isupp,work,iwork)

    !EE and BB plms and eigenvalues
    allocate(evec(nqu,nqu),eval(nqu),isupp(2*nqu),work(1),iwork(1))
    call dsyevr('V', 'A', 'L', nqu, Z(ntemp+1:ntot,1+ntemp:ntot), nqu, 0.d0, 0.d0, 0, 0, dlamch('S'), neigenl, eval, evec, nqu, isupp, work, -1, iwork, -1, info)
    lwork = work(1) ;liwork = iwork(1)
    deallocate(work,iwork) ;allocate(work(lwork),iwork(liwork))

    Z  = 0.d0
    do l =2,lswitch
       write(*,*) 'computing basis for PP l=',l
       il = nqu +1 -(2*l+1)
       iu = nqu
       
       clstmp = 0._dp
       clstmp(l,smwEE) = 1._dp
       call get_pp_cov(clstmp,l,l,Z(1+ntemp:ntot,1+ntemp:ntot),project_mondip=.false.)
       call dsyevr('V', 'I', 'L', nqu, Z(1+ntemp:ntot,1+ntemp:ntot), nqu, 0.d0, 0.d0, il, iu, 2*dlamch('S'), neigenl, eval, evec, nqu, isupp, work, lwork, iwork, liwork, info)
       eigenvec(l)%m(ntemp+1:ntot,neigenl+1:2*neigenl) = evec(1:nqu,1:neigenl)
       eigenval(l)%m(neigenl+1:2*neigenl,1)            = eval(1:neigenl)

       clstmp = 0._dp
       clstmp(l,smwBB) = 1._dp
       call get_pp_cov(clstmp,l,l,Z(1+ntemp:ntot,1+ntemp:ntot),project_mondip=.false.)
       call dsyevr('V', 'I', 'L', nqu, Z(1+ntemp:ntot,1+ntemp:ntot), nqu, 0.d0, 0.d0, il, iu, 2*dlamch('S'), neigenl, eval, evec, nqu, isupp, work, lwork, iwork, liwork, info)
       eigenvec(l)%m(ntemp+1:ntot,2*neigenl+1:3*neigenl) = evec(1:nqu,1:neigenl)
       eigenval(l)%m(2*neigenl+1:3*neigenl,1) = eval(1:neigenl)

    end do

    deallocate(evec,eval,isupp,work,iwork)

    !now find TE coefficients
    Z  = 0.d0
    do l=2,lswitch
       write(*,*) 'computing basis for TP l=',l

       clstmp = 0._dp
       clstmp(l,smwTE) = 1._dp
       call get_xx_cov(clstmp,l,l,Z,project_mondip=.false.,symmetrize=.true.)

       !allocate(evec,source=fevec(l)%m)
       allocate(evec(lbound(fevec(l)%m,1):ubound(fevec(l)%m,1), lbound(fevec(l)%m,2):ubound(fevec(l)%m,2)) &
        , source=fevec(l)%m)
       blocks(l)%m = 0._dp
       call dgemm('N','N', ntot, 2*(2*l+1), ntot, 1.d0, Z, ntot,&
            fevec(l)%m(:,1:2*(2*l+1)), ntot, 0.d0, evec ,ntot)
       call dgemm('T','N', 2*(2*l+1), 2*(2*l+1), ntot, 1.d0, fevec(l)%m(:,1:2*(2*l+1)), &
            ntot, evec, ntot, 0.d0, blocks(l)%m(1:2*(2*l+1),1:2*(2*l+1)) , 2*(2*l+1))

       deallocate(evec)
    end do

    deallocate(clstmp,Z)
    
    if (present(file)) call write_plms_to_file(file,eigenval,eigenvec,blocks)

  end subroutine compute_plms

  logical function read_plms_from_file(file,eigenval,eigenvec,blocks) result(OK)
    implicit none
    character(len=*) ,intent(in)            :: file
    type(allocatable_matrix) ,intent(inout) :: eigenval(2:) ,eigenvec(2:) &
         ,blocks(2:)
    integer(i4b)                            :: lsw ,flsw ,l ,fntot ,nl ,istat &
         ,unit

    unit = 3001
    lsw = size(eigenvec(:)) +1
    open(unit,file=trim(file),status='old',action='read',form='unformatted',iostat=istat)
    ok = istat .eq. 0
    if(.not. ok) return

    read(unit,iostat=istat) flsw
    ok = ok .and. (istat .eq. 0 .and. flsw .ge. lsw)

    read(unit,iostat=istat) fntot
    ok = ok .and. (istat .eq. 0 .and. fntot .eq. ntot)

    do l=2,lsw
       nl = 3*(2*l+1)
       read(unit,iostat=istat) eigenval(l)%m(1:nl,1)
       ok = ok .and. (istat .eq. 0)
       read(unit,iostat=istat) eigenvec(l)%m(1:fntot,1:nl)
       ok = ok .and. (istat .eq. 0)
       read(unit,iostat=istat) blocks(l)%m(1:nl,1:nl)
       ok = ok .and. (istat .eq. 0)
    end do
    close(unit)

  end function read_plms_from_file

  subroutine write_plms_to_file(file,eigenval,eigenvec,blocks)
    implicit none
    character(len=*) ,intent(in)         :: file
    type(allocatable_matrix) ,intent(in) :: eigenval(2:) ,eigenvec(2:) ,blocks(2:)
    integer(i4b)                         :: flsw ,l ,fntot ,nl ,unit

    unit = 3001

    open(unit,file=trim(file),status='unknown',action='write',form='unformatted')
    flsw  = size(eigenvec(:)) +1
    fntot = size(eigenvec(2)%m(:,1))
    write(unit) flsw
    write(unit) fntot
    do l=2,flsw
       nl = 3*(2*l+1)
       write(unit) eigenval(l)%m(1:nl,1)
       write(unit) eigenvec(l)%m(1:fntot,1:nl)
       write(unit) blocks(l)%m(1:nl,1:nl)
    end do
    close(unit)

  end subroutine write_plms_to_file

  subroutine clean_pix_like_smw
    implicit none

    integer :: l

    if(allocated(S))        deallocate(S)
    if(allocated(cls))      deallocate(cls)
    if(allocated(dt))       deallocate(dt)
    if(allocated(auxdt))    deallocate(auxdt)

    call ncov%clean()

    if(allocated(cblocks)) then
       do l=2,lswitch
          call cblocks(l)%clean()
       end do
       deallocate(cblocks)
    end if

    if(allocated(feval)) then
       do l=2,lswitch
          call feval(l)%clean()
       end do
       deallocate(feval)
    end if

  end subroutine clean_pix_like_smw

  subroutine get_tt_cov(clsin,linf,lsup,cov,project_mondip,symmetrize)

    implicit none
    real(dp) ,intent(in)          :: clsin(2:,:)
    integer ,intent(in)           :: linf ,lsup
    real(dp) ,intent(out)         :: cov(:,:)
    logical ,optional ,intent(in) :: project_mondip ,symmetrize

    integer  :: i ,j ,l
    real(dp) :: tt ,cz ,ct0 ,ct1

    ct0 = 0._dp
    ct1 = 0._dp

    if(present(project_mondip)) then
       if(project_mondip) then
          ct0 = 1.e6_dp
          ct1 = 1.e6_dp
       end if
    end if

    cls(2:lsup,smwTT) =clsin(2:lsup,smwTT)*clnorm(2:lsup,smwTT)

    do i=1,ntemp
! TT 
       cov(i,i) = sum(cls(linf:lsup,smwTT)) +ct1 +ct0        
       do j=i+1,ntemp
          cz = sum(Tvec(:,j)*Tvec(:,i))
          pl(1) = cz
          pl(2) = 1.5d0*cz*cz -.5d0
          do l = 3,lsup
             pl(l) =(cz*(2*l -1)*pl(l-1) -(l-1)*pl(l-2))/l
          enddo
          tt = sum(cls(linf:lsup,smwTT)*pl(linf:lsup)) &
               +Pl(1)*ct1 +ct0 !project out monopole and dipole
          cov(j,i) = tt
       enddo
    end do

    if (present(symmetrize) .and. symmetrize) then
       l = size(cov(:,1))
       do i=1,l
          do j=1,i-1
             cov(j,i) = cov(i,j)
          end do
       end do
    end if

  end subroutine get_tt_cov

  subroutine get_pp_cov(clsin,linf,lsup,cov,project_mondip,symmetrize)

    implicit none
    real(dp) ,intent(in) :: clsin(2:,:)
    integer ,intent(in)  :: linf ,lsup
    real(dp) ,intent(out):: cov(:,:)
    logical ,optional ,intent(in) :: project_mondip ,symmetrize

    integer  :: i,j,l,iu,iq,ju,jq
    real(dp) :: qq,uu,qu,cz
    real(dp) :: c1c2,s1s2,c1s2,s1c2

    cls(2:lsup,smwEE) = clsin(2:lsup,smwEE)*clnorm(2:lsup,smwEE)
    cls(2:lsup,smwBB) = clsin(2:lsup,smwBB)*clnorm(2:lsup,smwBB)
    cls(2:lsup,smwEB) = clsin(2:lsup,smwEB)*clnorm(2:lsup,smwEB)

    do i = 1,nq
!QQ 
       iq = i +ntemp
       do j = i,nq
          cz = sum(Qvec(:,j)*Qvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+cz*cz)
          f2(2)  = -12.d0*cz
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -cz*cz))*plm(l)+ &
                  (2*l+4)*cz*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*cz*plm(l) +(l+2)*plm(l-1))
          enddo
          qq = sum(cls(linf:lsup,smwEE)*f1(linf:lsup) -cls(linf:lsup,smwBB)*f2(linf:lsup))
          uu = sum(cls(linf:lsup,smwBB)*f1(linf:lsup) -cls(linf:lsup,smwEE)*f2(linf:lsup))
          qu = sum((f1(linf:lsup) +f2(linf:lsup))*cls(linf:lsup,smwEB))

          jq = j +ntemp

          c1c2 = cos1%column(iq)%row(jq)*cos2%column(iq)%row(jq)
          s1s2 = sin1%column(iq)%row(jq)*sin2%column(iq)%row(jq)
          c1s2 = cos1%column(iq)%row(jq)*sin2%column(iq)%row(jq)
          s1c2 = sin1%column(iq)%row(jq)*cos2%column(iq)%row(jq)

          cov(j,i) = qq*c1c2 +uu*s1s2 +qu*(c1s2+s1c2)

       end do

!UQ 
       do j=1,nu
          cz = sum(Uvec(:,j)*Qvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+cz*cz)
          f2(2)  = -12.d0*cz
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -cz*cz))*plm(l)+ &
                  (2*l+4)*cz*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*cz*plm(l) +(l+2)*plm(l-1))
          enddo
          qq = sum(cls(linf:lsup,smwEE)*f1(linf:lsup) -cls(linf:lsup,smwBB)*f2(linf:lsup))
          uu = sum(cls(linf:lsup,smwBB)*f1(linf:lsup) -cls(linf:lsup,smwEE)*f2(linf:lsup))
          qu = sum((f1(linf:lsup) +f2(linf:lsup))*cls(linf:lsup,smwEB))

          ju = j+nq+ntemp

          c1c2 = cos1%column(iq)%row(ju)*cos2%column(iq)%row(ju)
          s1s2 = sin1%column(iq)%row(ju)*sin2%column(iq)%row(ju)
          c1s2 = cos1%column(iq)%row(ju)*sin2%column(iq)%row(ju)
          s1c2 = sin1%column(iq)%row(ju)*cos2%column(iq)%row(ju)

          cov(j+nq,i) = -qq*s1c2 +uu*c1s2 +qu*(c1c2 -s1s2)
       end do
    end do

    do i =1,nu
!UU 
       iu = i+nq+ntemp
       do j=i,nu
          cz = sum(Uvec(:,j)*Uvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+cz*cz)
          f2(2)  = -12.d0*cz
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -cz*cz))*plm(l)+ &
                  (2*l+4)*cz*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*cz*plm(l) +(l+2)*plm(l-1))
          enddo
          qq = sum(cls(linf:lsup,smwEE)*f1(linf:lsup) -cls(linf:lsup,smwBB)*f2(linf:lsup))
          uu = sum(cls(linf:lsup,smwBB)*f1(linf:lsup) -cls(linf:lsup,smwEE)*f2(linf:lsup))
          qu = sum((f1(linf:lsup) +f2(linf:lsup))*cls(linf:lsup,smwEB))

          ju = j+nq+ntemp

          c1c2 = cos1%column(iu)%row(ju)*cos2%column(iu)%row(ju)
          s1s2 = sin1%column(iu)%row(ju)*sin2%column(iu)%row(ju)
          c1s2 = cos1%column(iu)%row(ju)*sin2%column(iu)%row(ju)
          s1c2 = sin1%column(iu)%row(ju)*cos2%column(iu)%row(ju)

          cov(j+nq,i+nq) = qq*s1s2 +uu*c1c2 -qu*(c1s2 +s1c2)
       end do
    end do

    if (present(symmetrize) .and. symmetrize) then
       l = size(cov(:,1))
       do i=1,l
          do j=1,i-1
             cov(j,i) = cov(i,j)
          end do
       end do
    end if

  end subroutine get_pp_cov

  subroutine get_xx_cov(clsin,linf,lsup,cov,project_mondip,symmetrize)

    implicit none
    real(dp) ,intent(in) :: clsin(2:,:)
    integer ,intent(in)  :: linf ,lsup
    real(dp) ,intent(out):: cov(:,:)
    logical ,optional ,intent(in) :: project_mondip ,symmetrize

    integer  :: i,j,l,ju,jq
    real(dp) :: tq,tu,cz

    cls(2:lsup,smwTE) =clsin(2:lsup,smwTE)*clnorm(2:lsup,smwTE)
    cls(2:lsup,smwTB) =clsin(2:lsup,smwTB)*clnorm(2:lsup,smwTB)


    do i=1,ntemp
! QT
       do j=1,nq
          cz = sum(Qvec(:,j)*Tvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0*(1.d0 -cz*cz)
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l-2)
          enddo
          tq = -sum(cls(linf:lsup,smwTE)*plm(linf:lsup))
          tu = -sum(cls(linf:lsup,smwTB)*plm(linf:lsup))

          jq = j +ntemp

          cov(jq,i) = tq*cos1%column(i)%row(jq)  +tu*sin1%column(i)%row(jq)

       enddo

!UT
       do j=1,nu
          cz = sum(Uvec(:,j)*Tvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0*(1.d0 -cz*cz)
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l-2)
          enddo
          tq = -sum(cls(linf:lsup,smwTE)*plm(linf:lsup))
          tu = -sum(cls(linf:lsup,smwTB)*plm(linf:lsup))

          ju = j +ntemp +nq

          cov(ju,i) = -tq*sin1%column(i)%row(ju) +tu*cos1%column(i)%row(ju)

       enddo

    end do

    if (present(symmetrize) .and. symmetrize) then
       l = size(cov(:,1))
       do i=1,l
          do j=1,i-1
             cov(j,i) = cov(i,j)
          end do
       end do
    end if

  end subroutine get_xx_cov

  subroutine get_pix_loglike_smw(clsin,alike)
!input clsin are l(l+1)C_l/2pi

    implicit none
    
    real(dp),intent(in)   :: clsin(2:,1:)
    real(dp),intent(out)  :: alike

    real(dp)              :: argexp,logdet

    integer               :: i ,j ,l ,info ,neigen ,nl ,nl2 ,nl3 ,iu
    integer ,allocatable  :: pinfo(:) ,il(:)
    real(dp) ,allocatable :: tm1(:,:) ,tmplog(:)
    real(dp)              :: tmp ,t1 ,t2

    neigen = size(dt(:,1))
    logdet = 0._dp

    nl2 = 2*(2*lswitch+1)
    allocate(tm1(nl2,nl2),tmplog(2:lswitch),pinfo(2:lswitch),il(2:lswitch))
    iu = 0
    do l=2,lswitch
       il(l) = iu +1
       iu = il(l) +6*l+2 !3*(2l+1) for each l
    end do

    call cpu_time(t1)
    do i=1,neigen
       do j=i,neigen
          S(j,i) = ncov%column(i)%row(j)
       end do
    end do

    do l=2,lswitch

       tmplog(l) = 0._dp
    
       nl = 2*l+1
       nl2 = 2*nl
       nl3 = 3*nl
 
       ! TT,TE,EE 
       do i=1,nl
          tm1(i,i)        = clsin(l,smwTT)*feval(l)%m(i,1)
          tm1(i+1:nl,i)   = 0.d0
          tm1(nl+1:nl2,i) = clsin(l,smwTE)*cblocks(l)%m(nl+1:nl2,i)
       end do
       do i = nl+1,nl2
          tm1(i,i)        = clsin(l,smwEE)*feval(l)%m(i,1)
          tm1(i+1:nl2,i)  = 0.d0
       end do

       call dpotrf('L',nl2,tm1(1:nl2,1:nl2),nl2,pinfo(l))

       do j=1,nl2
          tmplog(l) = tmplog(l) +log(tm1(j,j))
       end do
       call dpotri('L',nl2,tm1(1:nl2,1:nl2),nl2,info)

       do i=1,nl2
          do j=i,nl2
             S(il(l)-1+j,il(l)-1+i) = S(il(l)-1+j,il(l)-1+i) +tm1(j,i)
          end do
       end do

       !BB is diagonal 
       do i = nl2+1,nl3
          tmp = 1._dp/(max(clsin(l,smwBB),1.d-30)*feval(l)%m(i,1))
          S(il(l)+i-1,il(l)+i-1) = S(il(l)+i-1,il(l)+i-1) +tmp
          tmplog(l) = tmplog(l) -.5_dp*log(tmp) 
       end do
 
    end do
    call cpu_time(t2)

    !print*,'build matrix in :',t2-t1
    logdet = sum(tmplog(2:lswitch))
    
    do l=2,lswitch
       if(pinfo(l).ne.0) then
          write(*,*) 'block inversion failed for l =',l,' info =',pinfo(l)
          argexp = -1.d30
          logdet = -1.d30
          alike  = -1.d30
          deallocate(tmplog,il,tm1,pinfo)
          return
       end if
    end do

    call cpu_time(t1)
    auxdt = dt
    call dposv('L',neigen,ndata,S,neigen,auxdt,neigen,info)
    call cpu_time(t2)
    !print*,'dposv in :',t2-t1

!    cholesky: X*(C^-1*X)
    if (info.eq.0) then

       do j=1,neigen
          logdet = logdet +log(S(j,j))
       enddo
       logdet = 2.d0*logdet

       !remeber: this version for ndata =1
       argexp = reflike(1) -sum(dt(:,1)*auxdt(:,1))
       alike  = -.5d0*(argexp+logdet)

    else

       argexp = 1.d30
       logdet = 1.d30

       alike = -.5d0*(argexp+logdet)
    endif

    deallocate(tmplog,il,tm1,pinfo)
  end subroutine get_pix_loglike_smw

  subroutine update_ncvm(clsin,linf,lsup,NCM,project_mondip)
!input clsin are l(l+1)C_l/2pi

    implicit none
    
    real(dp)  ,intent(in)         :: clsin(2:,:)
    integer  ,intent(in)          :: linf ,lsup
    real(dp) ,intent(inout)       :: NCM(:,:)
    logical ,intent(in) ,optional :: project_mondip

    integer  :: i ,j ,l ,iu ,iq ,ju ,jq
    real(dp) :: tt ,qq ,uu ,tq ,tu ,qu ,cz ,c1c2 ,s1s2 ,c1s2 ,s1c2 ,ct0 ,ct1

    ct0 = 0._dp
    ct1 = 0._dp

    if(present(project_mondip)) then
       if(project_mondip) then
          ct0 = 1.e6_dp
          ct1 = 1.e6_dp
       end if
    end if

    cls(2:lsup,:) =clsin(2:lsup,:)*clnorm(2:lsup,:)

    do i=1,ntemp
! TT 
       do j=i,ntemp
          cz = sum(Tvec(:,j)*Tvec(:,i))
          pl(1) = cz
          pl(2) = 1.5d0*cz*cz -.5d0
          do l = 3,lsup
             pl(l) =(cz*(2*l -1)*pl(l-1) -(l-1)*pl(l-2))/l
          enddo
          tt = sum(cls(linf:lsup,smwTT)*pl(linf:lsup)) 
          NCM(j,i) = NCM(j,i) +tt +ct0 +pl(1)*ct1
       enddo

! QT
       do j=1,nq
          cz = sum(Qvec(:,j)*Tvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0*(1.d0 -cz*cz)
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l-2)
          enddo
          tq = -sum(cls(linf:lsup,smwTE)*plm(linf:lsup))
          tu = -sum(cls(linf:lsup,smwTB)*plm(linf:lsup))

          jq = j +ntemp

          NCM(jq,i) = NCM(jq,i) +tq*cos1%column(i)%row(jq) +tu*sin1%column(i)%row(jq)
       enddo

!UT
       do j=1,nu
          cz = sum(Uvec(:,j)*Tvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0*(1.d0 -cz*cz)
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l-2)
          enddo
          tq = -sum(cls(linf:lsup,smwTE)*plm(linf:lsup))
          tu = -sum(cls(linf:lsup,smwTB)*plm(linf:lsup))

          ju = j +ntemp +nq

          NCM(ju,i) = NCM(ju,i) -tq*sin1%column(i)%row(ju) +tu*cos1%column(i)%row(ju)
       enddo

    end do

    do i = 1,nq
!QQ 
       do j = i,nq
          cz = sum(Qvec(:,j)*Qvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+cz*cz)
          f2(2)  = -12.d0*cz
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -cz*cz))*plm(l)+ &
                  (2*l+4)*cz*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*cz*plm(l) +(l+2)*plm(l-1))
          enddo
          qq = sum(cls(linf:lsup,smwEE)*f1(linf:lsup) -cls(linf:lsup,smwBB)*f2(linf:lsup))
          uu = sum(cls(linf:lsup,smwBB)*f1(linf:lsup) -cls(linf:lsup,smwEE)*f2(linf:lsup))
          qu = sum((f1(linf:lsup) +f2(linf:lsup))*cls(linf:lsup,smwEB))

          jq = j+ntemp
          iq = i+ntemp

          c1c2 = cos1%column(iq)%row(jq)*cos2%column(iq)%row(jq)
          s1s2 = sin1%column(iq)%row(jq)*sin2%column(iq)%row(jq)
          c1s2 = cos1%column(iq)%row(jq)*sin2%column(iq)%row(jq)
          s1c2 = sin1%column(iq)%row(jq)*cos2%column(iq)%row(jq)

          NCM(jq,iq) = NCM(jq,iq) +qq*c1c2 +uu*s1s2 +qu*(c1s2+s1c2)

       end do

!UQ 
       do j=1,nu

          cz = sum(Uvec(:,j)*Qvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+cz*cz)
          f2(2)  = -12.d0*cz

          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -cz*cz))*plm(l)+ &
                  (2*l+4)*cz*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*cz*plm(l) +(l+2)*plm(l-1))
          enddo

          qq = sum(cls(linf:lsup,smwEE)*f1(linf:lsup) -cls(linf:lsup,smwBB)*f2(linf:lsup))
          uu = sum(cls(linf:lsup,smwBB)*f1(linf:lsup) -cls(linf:lsup,smwEE)*f2(linf:lsup))
          qu = sum((f1(linf:lsup) +f2(linf:lsup))*cls(linf:lsup,smwEB))

          ju = j+ntemp+nq
          iq = i+ntemp

          c1c2 = cos1%column(iq)%row(ju)*cos2%column(iq)%row(ju)
          s1s2 = sin1%column(iq)%row(ju)*sin2%column(iq)%row(ju)
          c1s2 = cos1%column(iq)%row(ju)*sin2%column(iq)%row(ju)
          s1c2 = sin1%column(iq)%row(ju)*cos2%column(iq)%row(ju)

          NCM(ju,iq) = NCM(ju,iq) -qq*s1c2 +uu*c1s2 +qu*(c1c2 -s1s2)
       end do
    end do

    do i =1,nu
!UU 
       do j=i,nu
          cz = sum(Uvec(:,j)*Uvec(:,i))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+cz*cz)
          f2(2)  = -12.d0*cz
          do l = 3,lsup
             plm(l) =(cz*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -cz*cz))*plm(l)+ &
                  (2*l+4)*cz*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*cz*plm(l) +(l+2)*plm(l-1))
          enddo
          qq = sum(cls(linf:lsup,smwEE)*f1(linf:lsup) -cls(linf:lsup,smwBB)*f2(linf:lsup))
          uu = sum(cls(linf:lsup,smwBB)*f1(linf:lsup) -cls(linf:lsup,smwEE)*f2(linf:lsup))
          qu = sum((f1(linf:lsup) +f2(linf:lsup))*cls(linf:lsup,smwEB))


          ju = j+ntemp+nq
          iu = i+ntemp+nq

          c1c2 = cos1%column(iu)%row(ju)*cos2%column(iu)%row(ju)
          s1s2 = sin1%column(iu)%row(ju)*sin2%column(iu)%row(ju)
          c1s2 = cos1%column(iu)%row(ju)*sin2%column(iu)%row(ju)
          s1c2 = sin1%column(iu)%row(ju)*cos2%column(iu)%row(ju)

          NCM(ju,iu) = NCM(ju,iu) +qq*s1s2 +uu*c1c2 -qu*(c1s2 +s1c2)
       end do
    end do

  end subroutine update_ncvm

  pure subroutine get_rotation_angle(r1,r2,a12,a21)
!computes TWICE the rotation angle

    implicit none
    real(dp) ,intent(in) ,dimension(3) :: r1,r2
    real(dp) ,intent(out)              :: a12
    real(dp) ,intent(out) ,optional    :: a21

    real(dp) ,parameter ::eps = 3.141592653589793d0/180.d0/3600.d0/100.d0
    real(dp) ,parameter ,dimension(3) :: zz =(/0,0,1/),epsilon =(/eps,0.d0,0.d0/)
    
    real(dp) ,dimension(3) :: r12,r1star,r2star
    real(dp)               :: mod

    call ext_prod(r1,r2,r12)
    mod = sqrt(sum(r12*r12))
    if(mod.lt.1.d-8) then !same or opposite pixels
       a12 = 0.d0
       if(present(a21)) a21 = 0.d0
       return
       
    end if
    r12 = r12/mod

    call ext_prod(zz,r1,r1star)
    r1star(3) = 0.d0
    mod = sqrt(sum(r1star*r1star))
    if(mod.lt.1.d-8) then   !r1 is at a pole            
       r1star = r1star+epsilon
       mod = sqrt(sum(r1star*r1star))
    end if
    r1star = r1star/mod

    call ext_prod(zz,r2,r2star)
    r2star(3) = 0.d0
    mod = sqrt(sum(r2star*r2star))
    if(mod.lt.1.d-8) then   !r2 is at a pole            
       r2star = r2star+epsilon
       mod = sqrt(sum(r2star*r2star))
    end if
    r2star = r2star/mod

    mod = sum(r12*r1star)
    mod = min(1.d0,mod)
    mod = max(-1.d0,mod)
    if(sum(r12*zz).gt.0.d0) then
       a12 = 2.d0*acos(mod)
    else
       a12 = -2.d0*acos(mod)
    end if

    a12 = -a12

    if(present(a21)) then
       r12 = -r12   !r21 = - r12
       mod = sum(r12*r2star)
       mod = min(1.d0,mod)
       mod = max(-1.d0,mod)
       if(sum(r12*zz).gt.0.d0) then
          a21 =  2.d0*acos(mod)
       else
          a21 = -2.d0*acos(mod)
       end if

       a21 = -a21
    end if

  end subroutine get_rotation_angle

  pure subroutine ext_prod(v1,v2,v3)
    implicit none
    real(dp),intent(in) :: v1(3),v2(3)
    real(dp),intent(out) :: v3(3)
    
    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

  end subroutine ext_prod

  logical function read_cls_from_file(filein,cls) result(ok)

    use fiducial_cls_mod_smw
    implicit none

    character(len=*) ,intent(in) :: filein
    real(dp) ,intent(out)        :: cls(2:,:)

    integer   :: istat ,lmax ,unit ,l ,idum

    lmax = size(cls(:,1)) +1

    cls = 0._dp
    ok = .false.

    call fid%init()

    open(unit = unit,file=trim(filein),status='old',action='read',&
         iostat=istat)
    if(istat .ne. 0) then
       cls(2:lmax,6) = fid%cls(2:lmax,6)
       return
    end if
    read(unit,*,iostat=istat) idum,cls(2,1:4)
    if(istat .eq. 0) then 
       write(*,*) 'cls file appears to have 5+ columns'
       write(*,*) 'assuming it is a CAMB file with l, TT, EE, BB, TE'
       do l=3,lmax
          read(unit,*,iostat=istat) idum,cls(l,1:4)
       end do
    else
       rewind(unit)
       read(unit,*,iostat=istat) idum,cls(2,smwTT),cls(2,smwEE),cls(2,smwTE)
       if(istat .eq.0) then
          write(*,*) 'cls file appears to have 4 columns'
          write(*,*) 'assuming it is a CAMB file with l, TT, EE, TE'
          do l=3,lmax
             read(unit,*,iostat=istat) idum,cls(l,smwTT),cls(l,smwEE),cls(l,smwTE)
          end do
       else
          cls(2:lmax,6) = fid%cls(2:lmax,6)
          return
       end if
    end if
    ok = .true.

  end function read_cls_from_file

end module 



