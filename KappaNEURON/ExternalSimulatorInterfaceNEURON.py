import os
import re
import SpatialKappa
from py4j.protocol import Py4JJavaError
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
        self._species = kwargs.get('species', [])
        kappa_file = kwargs.get('kappa_file')
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
            None,  # ode_count
            None,  # ode_reinit
            None,  # rhs
            None,  # ode_solve
            None,  # ode_jacobian
            None,  # ode_abs_tolerance
        ]
        self.register()

        nodelist = self._species[0][self._region].nodes
        self._segments = [n.segment for n in nodelist]

    def init(self):
        global gateway

        self._kappa_sims = []  # Will this destroy things properly?
        self._name = dict()

        for segment in self._segments:
            kappa_sim = gateway.kappa_sim(self._time_units, True, self.seed)
            self._name[segment] = kappa_sim

            try:
                kappa_sim.loadFile(self._kappa_file)
            except Py4JJavaError as e:
                java_err = re.sub(r'java.lang.IllegalStateException: ', r'', str(e.java_exception))
                errstr = 'Error in kappa file %s: %s' % (self._kappa_file, java_err)
                raise RuntimeError(errstr)
            kappa_sim.initialiseSim()
            self._kappa_sims.append(kappa_sim)

    def run_every_fixed_timestep(self, dt):

        for kappa_sim in self._kappa_sims:
            kappa_sim.runForTime(dt, False)

        for s in self._species:
            nodelist = s[self._region].nodes
            for n in nodelist:
                kappa_sim = self._name[n.segment]
                n.concentration = kappa_sim.getObservation(s.name) #/ (n.volume * 6.022 * pow(10, 5))

    def register(self):
        nbs.register(self._callbacks)

    def unregister(self):
        nbs.unregister(self._callbacks)


dend = h.Section(name="dend")
dend.nseg = 3
cyt = rxd.Region([dend], name="cyt")
A = rxd.Species([cyt], name="A", charge=0, initial=0 * uM)
B = rxd.Species([cyt], name="B", charge=0, initial=0 * uM)
AB = rxd.Species([cyt], name="AB", charge=0, initial=0 * uM)

t = h.Vector().record(h._ref_t)
A_cyt01 = h.Vector().record(A[cyt].nodes(dend(0.1))._ref_concentration)
B_cyt01 = h.Vector().record(B[cyt].nodes(dend(0.1))._ref_concentration)
AB_cyt01 = h.Vector().record(AB[cyt].nodes(dend(0.1))._ref_concentration)

es1 = ExternalSimulatorInterface(cyt, species=[A, B, AB], kappa_file="ab.ka")

h.finitialize(-65 * mV)
h.continuerun(100 * ms)

A_cyt01[0] = 150
B_cyt01[0] = 100

plt.plot(t, A_cyt01, label="A")
plt.plot(t, B_cyt01, label="B")
plt.plot(t, AB_cyt01, label="AB")
plt.xlabel("Time(ms)")
plt.ylabel("Number of molecule")

plt.legend()
plt.show()
