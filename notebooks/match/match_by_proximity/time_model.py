from numcosmo_py import Ncm, GObject, Nc

"""Example of defining a new Ncm.Model for SLine."""

from numcosmo_py import Ncm, GObject

#
# New ModelBuilder object, defines a new model NcHIPrimExample
# implementing the base Ncm.Model abstract class.
#
mb = Ncm.ModelBuilder.new(Ncm.Model, "NcPySLineModel", "A simple python example model")

#
# New parameter `alpha' to describe the ln-slope
# Allowed interval: [0, 5]; Default scale: 0.1
# Absolute tolerance: 0; Default value: 2
#
mb.add_sparam(r"\alpha", "alpha", -10.0, 10.0, 0.1, 0.0, 2.0, Ncm.ParamType.FREE)

#
# New parameter `a' to describe the amplitude
# Allowed interval: [0.2, 2]; Default scale: 0.1
# Absolute tolerance: 0; Default value: 1
#
mb.add_sparam("a", "a", -5.0, 5.0, 0.1, 0.0, 1.0, Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet! Then register
# the new object in the GObject type system by creating a new
# instance. Finally, gets the Python version of the object (.pytype)
# and register it as a PyGObject object.
#
GNcPySLineModel = mb.create()
GObject.new(GNcPySLineModel)
NcPySLineModel = GNcPySLineModel.pytype
GObject.type_register(NcPySLineModel)

class PySLineModel(NcPySLineModel):  # type: ignore
    """A simple python example model."""

    #
    # Defining some property which is not part of the model paramers.
    # All model parameters must be defined by the ModelBuilder.
    #
    some_property = GObject.Property(type=str)

    #
    # Calling the father's constructor
    #
    def __init__(self):
        NcPySLineModel.__init__(self)

    #
    # Method to calculate the y(x)
    #
    def f_x(self, x):
        """Method to calculate the y(x)."""
        return self.props.alpha * x + self.props.a




class PySLineGauss(Ncm.DataGaussCov):
    """A simple data class for the SLine model based on Ncm.DataGaussCov."""

    #
    # We need one vector property to save the independent variables x
    #
    xv = GObject.Property(type=Ncm.Vector, flags=GObject.PARAM_READWRITE)

    #
    # The contructor assigns some default values and calls the father's constructor.
    #
    def __init__(self, length=600):
        Ncm.DataGaussCov.__init__(self, n_points=length)
        n_points = super().get_size()

        if n_points > 0:
            self.xv = Ncm.Vector.new(n_points)
        else:
            self.xv = None

        self.cov_init = False

        #
        # Initializing to sane values
        #
        if n_points > 0:
            self.peek_cov().set_identity()
            self.xv.set_zero()

    #
    # Implements the virtual method get_length.
    #
    def do_get_length(self) -> int:  # pylint: disable-msg=arguments-differ
        return super().get_size()

    #
    # Implements the virtual method get_dof.
    #
    def do_get_dof(self) -> int:  # pylint: disable-msg=arguments-differ
        return super().get_size()

    #
    # Implements the virtual method `begin'.
    # This method usually do some groundwork in the data
    # before the actual calculations. For example, if the likelihood
    # involves the decomposition of a constant matrix, it can be done
    # during `begin' once and then used afterwards.
    #
    def do_begin(self):  # pylint: disable-msg=arguments-differ
        return

    #
    # Implements the virtual method `prepare'.
    # This method should do all the necessary calculations using mset
    # to be able to calculate the likelihood afterwards.
    #
    def do_prepare(self, _mset):  # pylint: disable-msg=arguments-differ
        return

    #
    # Implements the virtual method `mean_func'.
    # This method should compute the theoretical mean for the gaussian
    # distribution.
    #
    # pylint: disable-next=arguments-differ
    def do_mean_func(self, mset: Ncm.MSet, vp: Ncm.Vector) -> None:
        mid = mset.get_id_by_ns("NcPySLineModel")
        slm = mset.peek(mid)
        assert isinstance(slm, PySLineModel)
        n_points = super().get_size()

        for i in range(n_points):
            x = self.xv.get(i)
            vp.set(i, slm.f_x(x))




#
# Register our new Python class PySLineGauss
#
GObject.type_register(PySLineGauss)

#
# Register our new Python class PyNcPySLineModel
#
GObject.type_register(PySLineModel)


