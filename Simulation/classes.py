import math
import parameters
import numpy as np


class Coil:
    def __init__(self, material, r_inner, length, layers, wire_d):
        self.current = 0
        self.energy = 0
        self.temp = parameters.ambient_temp
        self.time_on = 0

        self.r_inner = r_inner
        self.length = length

        self.layers = layers
        self.wire_d = wire_d
        self.turns = round(length * layers / wire_d)
        self.r_outer = r_inner + (layers - 1) * wire_d
        self.r_avg = (r_inner + self.r_outer) / 2
        self.area = pow(self.r_avg, 2) * np.pi

        self.material = material
        self.wire_length = self.turns * np.pi * self.r_avg * 2
        self.wire_area = math.pi * pow(self.wire_d / 2, 2)
        self.wire_mass = self.wire_area * self.wire_length * material.density
        self.thermal_mass = self.wire_mass * material.heat_capacity
        self.wire_resistance = material.resistivity * self.wire_length / self.wire_area

        self.inductance = parameters.mu0 * self.area * pow(self.turns, 2) / length
        self.resistance = self.wire_resistance
        self.tau = self.inductance / self.resistance

        self.currents_list = []
        self.energy_list = []

        print("Wire Length", self.wire_length, "m")
        print("Wire Mass", self.wire_mass, "kg")
        print("Turns", self.turns)
        print("mass", self.wire_mass, "kg")


class Stage:
    def __init__(self, coil, position, trigger_position, capacitance, voltage):
        self.coil = coil
        self.pos = position  # Refers to axial distance of coil center
        self.trigger_pos = trigger_position

        self.capacitance = capacitance
        self.voltage = voltage

        coil.resistance = 2 * math.sqrt(coil.inductance / self.capacitance)
        coil.added_resistance = coil.resistance - coil.wire_resistance
        coil.pos = position
        print("Added resistance", coil.added_resistance, "Ohms")


class Projectile:
    def __init__(self, coil, position, mass):
        self.coil = coil
        self.pos = position  # Refers to axial distance of coil center
        self.vel = 0
        self.kinetic_energy = 0
        self.mass = mass
        self.stage = 0

        self.position_list = []
        self.velocity_list = []
        self.acceleration_list = []
        self.energy_list = []


class Material:
    def __init__(self, density, heat_capacity, resistivity, alpha):
        self.density = density
        self.heat_capacity = heat_capacity
        self.resistivity = resistivity
        self.alpha = alpha

