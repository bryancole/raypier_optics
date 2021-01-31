
import unittest

from raypier.tracer import RayTraceModel
from raypier.waveplates import Waveplate
from raypier.sources import ParallelRaySource, SingleRaySource
from raypier.mirrors import PECMirror

from numpy import pi
import numpy


class TestWaveplate(unittest.TestCase):
    def setUp(self):
        self.wp = Waveplate(direction=(0.,0.,1.))
        self.src = ParallelRaySource(rings=0,
                                     E_vector=(1,0,0),
                                     E1_amp = 1.0+0.0j,
                                     E2_amp = 0.0+0.0j,
                                     direction = (0,0,1),
                                     origin = (0,0,-10),
                                    )

        self.model = RayTraceModel(optics=[self.wp,],
                                   sources=[self.src])

    def test_retardance_property(self):
        wp = self.wp
        for val in [0.0, 0.25, 0.5, 0.75, 1.0]:
            wp.retardance = val
            self.assertEqual(wp.retardance, val)

    def test_input_ray(self):
        ray = self.src.input_rays[0]
        print(ray.power, ray.ellipticity, ray.E_vector, ray.E1_amp, ray.E2_amp)
        self.assertEqual(ray.ellipticity, 0.0)

    def test_zero_retardance(self):
        wp = self.wp
        wp.retardance = 0.0
        wp.rotation = 19.235
        direction=(0.,0.,1.)
        print(wp.x_axis)
        self.model.trace_all()

        ray_in = self.src.input_rays[0]
        ray_out = self.src.traced_rays[-1][0]

        ray_in.project_E(1,1,0)
        ray_out.project_E(1,1,0)

        self.assertEqual(ray_in.E_vector, ray_out.E_vector)
        self.assertAlmostEqual(ray_in.E1_amp, ray_out.E1_amp)
        self.assertAlmostEqual(ray_in.E2_amp, ray_out.E2_amp)

    def test_zero_retardance2(self):
        wp = self.wp
        wp.retardance = 0.0
        wp.rotation = 19.235
        direction=(0.,0.,-1.)
        print(wp.x_axis)
        self.model.trace_all()

        ray_in = self.src.input_rays[0]
        ray_out = self.src.traced_rays[-1][0]

        ray_in.project_E(1,1,0)
        ray_out.project_E(1,1,0)

        self.assertEqual(ray_in.E_vector, ray_out.E_vector)
        self.assertAlmostEqual(ray_in.E1_amp, ray_out.E1_amp)
        self.assertAlmostEqual(ray_in.E2_amp, ray_out.E2_amp)

    def test_half_wave_plate(self):
        wp = self.wp
        wp.retardance = 0.5
        for rot in numpy.linspace(-180,180,20):
            wp.rotation = rot
            self.model.trace_all()

            ray_in = self.src.input_rays[0]
            ray_out = self.src.traced_rays[-1][0]

            self.assertAlmostEqual(ray_in.ellipticity, ray_out.ellipticity)


class TestWaveplateBackreflection(unittest.TestCase):
    def test_quarter_retardance(self):
        src = SingleRaySource(E_vector=(1,0,0),
                                 E1_amp = 1.0+0.0j,
                                 E2_amp = 0.0+0.0j,
                                 direction = (0,0,1),
                                 origin = (0,0,-40),
                                )

        wp = Waveplate(retardance=0.25,
                       rotation = 45.0,
                       centre=(0,0,40),
                       direction=(0,0,1))

        m = PECMirror(centre=(0,0,80),
                      direction=(0.0,0.0,1.0))

        model = RayTraceModel(optics=[wp,m], sources=[src])
        #model.trace_all()

        def get_major_axis(ray):
            E1_vector = numpy.array(ray.E_vector)
            E2_vector = numpy.cross(numpy.array(ray.direction),E1_vector
                                    )
            E = (E1_vector*ray.E1_amp + E2_vector*ray.E2_amp)
            return E

        ray_in = src.input_rays[0]

        #print ray_in.E_vector, ray_in.E1_amp, ray_in.E2_amp

        self.assertEqual(len(src.traced_rays), 4)

        for i in range(4):
            ray = src.traced_rays[i][0]
            print(i, "::", ray.E_left, ray.E_right, ray.E_vector, ray.E1_amp, ray.E2_amp)
            print("\t", get_major_axis(ray))

        ray_out = src.traced_rays[-1][0]
        self.assertAlmostEqual(numpy.dot(get_major_axis(ray_out),
                                    get_major_axis(ray_in)), 0.0)

    def test_quarter_retardance2(self):
        src = SingleRaySource(E_vector=(1,0,0),
                                 E1_amp = 1.0+0.0j,
                                 E2_amp = 0.0+0.0j,
                                 direction = (0,0,1),
                                 origin = (0,0,-40),
                                )

        wp1 = Waveplate(retardance=0.25,
                       rotation = 45.0,
                       centre=(0,0,40),
                       direction=(0,0,1))

        print("wp1 fa:", wp1.material.fast_axis)

        wp2 = Waveplate(retardance=0.25,
                       rotation = -45.0,
                       centre=(0,0,80),
                       direction=(0,0,-1))

        print("wp2 fa:", wp2.material.fast_axis)
        print("wp1 fa:", wp1.material.fast_axis)

        model = RayTraceModel(optics=[wp1,wp2], sources=[src])
        model.trace_all()

        def get_major_axis(ray):
            E1_vector = numpy.array(ray.E_vector)
            E2_vector = numpy.cross(numpy.array(ray.direction),
                                    E1_vector)
            E = (E1_vector*ray.E1_amp + E2_vector*ray.E2_amp)
            return E.real

        ray_in = src.input_rays[0]

        #print ray_in.E_vector, ray_in.E1_amp, ray_in.E2_amp

        self.assertEqual(len(src.traced_rays), 3)

        for i in range(len(src.traced_rays)):
            ray = src.traced_rays[i][0]
            print(i, "::", get_major_axis(ray), ray.E_left, ray.E_right, ray.ellipticity, ray.E_vector)

if __name__ == "__main__":
    unittest.main()
