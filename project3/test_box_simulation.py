from Simulation_box import simulation_box 
import math 
class TestSimulationBox: 
    def __init__(self):
        self.box=simulation_box(0,20,0,20,0,20)

    def test_walker_simulation_no_sphere(self,steps):
        """ 
        Unit Test for the walker method, here we test on a simulation with no sphere, wich mean that the estimated volume of sphere should be around 0.
        The result could be slighly over 0 if the walkers are not going on all the map, that's why we add abs_error of 200m^3
        """
        print("\n --- Start of test with no spheres---\n")
        self.box = simulation_box(0,20,0,20,0,20)
        self.box.generate_sphere(n=0,specific=False)
        self.box.generate_starting_walker(n=2000)
        Monte_Carlo_volume=self.box.start_walker_simulation(steps=steps)
        assert math.isclose(Monte_Carlo_volume,0,abs_tol=200),"Volume should be equal to the box volume when no spheres are present"

    def test_walker_simulation_with_sphere(self,steps):
        """ 
        Unit Test for the walker method, here we test on a simulation with many spheres, wich mean that the estimated volume of sphere should be around the volume of the box.
        The result could be slighly less then 8000 because the distribution of sphere is random, that's why we add rel_error of 5% of the volume of the box.
        """
        print("\n--- Start of test with many spheres---\n")
        self.box = simulation_box(0,20,0,20,0,20)
        self.box.generate_sphere(n=200,specific=False)
        self.box.generate_starting_walker(n=2000)
        Monte_Carlo_volume=self.box.start_walker_simulation(steps=steps)
        print("Walker method estimation :",Monte_Carlo_volume)
        assert math.isclose(Monte_Carlo_volume,8000,rel_tol=0.05),"Volume Should be close to 8000 with that many spheres"

    def test_walker_simulation_random(self,steps,n_s=5):
        """ 
        Unit Test for the walker method, here we test on a simulation with a specific number of sphere.
        We calculate the result with walker method and we compare it to the analytical calculation from the data of the sphere. 
        The result might be different with big values of sphere cause the analytical calculation add each volume of sphere, 
        it count twice when they are on the same coordonate.
        This test is ne very stable
        """
        print("\n--- Start of test with random spheres---\n")
        self.box = simulation_box(0,20,0,20,0,20)
        self.box.generate_sphere(n=n_s,specific=False)
        self.box.generate_starting_walker(n=2000)
        walker_volume=self.box.start_walker_simulation(steps=steps)
        analytical_volume=self.box.test_volume_sphere()
        print("Estimated Volume :", walker_volume, "Analytical Volume :", analytical_volume)
        assert math.isclose(walker_volume,analytical_volume,rel_tol=0.15),"Walker Simulation is too far from the reality"
    def test_pi_approximation(self):
        """ 
        Unit test to try the estimation of PI method.
        """
        print("\n--- Start of test calculation of pi---\n")
        self.box = simulation_box(0,1,0,1,0,1)
        self.box.generate_sphere(n=1,specific=True)
        self.box.random_point(n=2000)
        self.box.check_point_spot()
        pi_estimate=self.box.estimate_pi()
        assert(pi_estimate>3 and pi_estimate<3.3), "The estimated value of pi is out of range !"

    def test_part_1_task_7_volume(self):
        """
        Unit test for the calculation of sphere volume, for the task 7
        """
        print("\n --- Start of test on Task 7  volume estimation---\n")
        self.box = simulation_box(0,20,0,20,0,20)
        self.box.generate_sphere(n=10,specific=False)
        self.box.random_point(n=8000)
        monte_carlo_volume=self.box.check_point_spot()
        analytical_volume=self.box.test_volume_sphere()
        print("Analytical solution :",analytical_volume,"Monte Carlo value :",monte_carlo_volume)
        assert math.isclose(analytical_volume,monte_carlo_volume,rel_tol=0.30),"Monte Carlo value is to far from the reality"

    def test_part_1_task_4_volume(self):
        """
        Unit test for the calculation of sphere volume, for the task 4
        """
        print("\n --- Start of test on Task 4 volume estimation---\n")
        self.box = simulation_box(0,20,0,20,0,20)
        self.box.generate_sphere(n=1,specific=False)
        self.box.random_point(n=5000)
        monte_carlo_volume=self.box.check_point_spot()
        analytical_volume=self.box.test_volume_sphere()
        print("Analytical solution :",analytical_volume,"Monte Carlo value :",monte_carlo_volume)
        assert math.isclose(analytical_volume,monte_carlo_volume,abs_tol=200),"Monte Carlo value is to far from the reality"

    def test_part_1_task_10_volume(self):
        """
        Unit test for the calculation of sphere volume, for the task 10
        """
        print("\n --- Start of test on Task 10 volume estimation---\n")
        self.box = simulation_box(-50,-30,-20,10,-10,10)
        self.box.read_dna_coords()
        self.box.generate_sphere()
        self.box.random_point(n=5000)
        monte_carlo_volume=self.box.check_point_spot()
        analytical_volume=self.box.test_volume_sphere()
        print("Analytical solution :",analytical_volume,"Monte Carlo value :",monte_carlo_volume)
        assert math.isclose(analytical_volume,monte_carlo_volume,rel_tol=0.20),"Monte Carlo value is to far from the reality"

def launch_test():
    test_box=TestSimulationBox()
    test_box.test_pi_approximation()
    test_box.test_part_1_task_4_volume()
    test_box.test_part_1_task_7_volume()
    test_box.test_part_1_task_10_volume()

    print("All Part1 tests succes")

def launch_walker_test():
    test_box=TestSimulationBox()
    test_box.test_walker_simulation_no_sphere(steps=1000)
    test_box.test_walker_simulation_with_sphere(steps=100)
    test_box.test_walker_simulation_random(steps=1000)
    print("All Walker tests succes")
