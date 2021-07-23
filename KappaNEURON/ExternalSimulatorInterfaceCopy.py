import os
import re
import SpatialKappa
from py4j.protocol import Py4JJavaError
from neuron.rxd.generalizedReaction import GeneralizedReaction
from neuron import h, rxd
from neuron import nonvint_block_supervisor as nbs
from neuron.units import mV, ms, uM
import matplotlib.pyplot as plt

h.load_file("stdrun.hoc")

gateway = SpatialKappa.SpatialKappa(redirect_stdout=None)

def report(mess):
    global verbose
    if (verbose):
        print(mess)

class ExternalSimulatorInterface():
    def __init__(self, region, *args, **kwargs):

        # additional keyword arguments
        # self.membrane_species = kwargs.get('membrane_species', [])
        # self.species = kwargs.get('species', [])
        kappa_file = kwargs.get('kappa_file')
        # self.regions = kwargs.get('regions', None)
        # self.membrane_flux = kwargs.get('membrane_flux', True)
        time_units = kwargs.get('time_units', 'ms')
        self.seed = kwargs.get('seed', None)
        self._sk_redirect_stdout = kwargs.get('sk_redirect_stdout', None)
        self._time_units = time_units
        self._kappa_file = os.path.join(os.getcwd(), kappa_file)
        self.kappa_sim = None

        self._region = region
        self._callbacks = [
            None,  # setup
            self.init,  # initialize
            None,  # current
            None,  # conductance
            self.run_every_fixed_timestep,  # fixed step advance
            None,  # ode_count (need for variable step)
            None,  # ode_reinit (need for variable step)
            None,  # rhs (need for variable step)
            None,  # ode_solve (actually invert the Jacobian; need for variable step)
            None,
            # ode_jacobian (need for variable step... this gets called when the old Jacobian was no longer giving good estimates)
            None,  # ode_abs_tolerance
        ]
        self.register()

    def init(self):
        # # in principle, pass along the command to initialize
        # print(f"The external simulator interface for node at {self._node.segment}"
        #       f" in {[r.name for r in self._node.region]} has been initialized.")
        global gateway
        # global verbose

        ## Start Java and load the SpatialKappa class, if not already
        ## loaded
        # if not gateway:
        #     gateway = SpatialKappa.SpatialKappa(redirect_stdout=self._sk_redirect_stdout)

        # self._kappa_sims = []  # Will this destroy things properly?
        # for index in self._indices_dict[self._involved_species[0]()]:
        #     report("Creating Kappa Simulation in index %d" % (index))
        self.kappa_sim = gateway.kappa_sim(self._time_units, True, self.seed)
        try:
            self.kappa_sim.loadFile(self._kappa_file)
        except Py4JJavaError as e:
            java_err = re.sub(r'java.lang.IllegalStateException: ', r'', str(e.java_exception))
            errstr = 'Error in kappa file %s: %s' % (self._kappa_file, java_err)
            raise RuntimeError(errstr)
        print("initialing Spatialkappa")
        self.kappa_sim.initialiseSim()

    def run_every_fixed_timestep(self, dt):
        print("run every....")
        # do the calculation for the given advance
        # here is where an external simulator might be called
        result = h.sin(h.t)
        # now transfer the value to the node(s)
        # NOTE: if your external tool is working in molecule counts, be sure to
        # .      convert; can always ask a node for its .volume
        # self._node.concentration = result

        self.kappa_sim.runForTime(dt, False)
        print(self.kappa_sim.getTime())


    def register(self):
        nbs.register(self._callbacks)

    def unregister(self):
        # unfortunately cannot just have this in a __del__ because the nbs keeps
        # track of the callbacks which implicitly keeps a reference to this object
        nbs.unregister(self._callbacks)

dend = h.Section(name="dend")
dend.nseg = 3
# er = rxd.Region([dend], name="er", geometry=rxd.FractionalVolume(0.17))
cyt = rxd.Region([dend], name="cyt")
# ca = rxd.Species([cyt, er], name="ca", charge=2, initial=0 * uM)
A = rxd.Species([cyt], name = "A", charge = 0, initial= 1 * uM )
B = rxd.Species([cyt], name = "B", charge = 0, initial= 1 * uM )
AB = rxd.Species([cyt], name = "AB", charge = 0, initial= 1 * uM )

t = h.Vector().record(h._ref_t)
A_cyt01 = h.Vector().record(A[cyt].nodes(dend(0.1))._ref_concentration)
B_cyt01 = h.Vector().record(B[cyt].nodes(dend(0.1))._ref_concentration)
AB_cyt01 = h.Vector().record(AB[cyt].nodes(dend(0.1))._ref_concentration)

# ca_cyt09 = h.Vector().record(ca[cyt].nodes(dend(0.9))._ref_concentration)
# ca_er01 = h.Vector().record(ca[er].nodes(dend(0.1))._ref_concentration)
# ca_er09 = h.Vector().record(ca[er].nodes(dend(0.9))._ref_concentration)

es1 = ExternalSimulatorInterface(cyt, kappa_file = "ab.ka")
# es2 = ExternalSimulatorInterface(ca[er].nodes(dend(0.9)), kappa_file = "tests/caPump1.ka")

h.finitialize(-65 * mV)
h.continuerun(1.5 * ms)

# es2.unregister()

h.continuerun(3 * ms)

subplots = []
fig = plt.figure(figsize=(12, 4))
for vec in [A_cyt01, B_cyt01, AB_cyt01]:
  subplot = fig.add_subplot(2, 2, len(subplots) + 1)
  subplot.plot(t, vec)
  subplot.set_ylim(-0.5, 1.5)
  # if vec in [ca_er01, ca_er09]:
  #   subplot.set_yticks([])
  subplots.append(subplot)

# subplots[0].set_ylabel("dend(0.1) conc")
# subplots[0].set_title("cyt")
# subplots[1].set_title("er")
# subplots[2].set_ylabel("dend(0.9) conc")
# subplots[2].set_xlabel("t (ms)")
# subplots[3].set_xlabel("t (ms)")

fig.show()
