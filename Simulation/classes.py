import math
import parameters
import numpy as np


class Coil:
    def __init__(self, material, r_inner, length, layers, wire_d):
        self.current = 0
        self.voltage = 0
        self.energy = 0
        self.temp = parameters.ambient_temp
        self.time_on = 0
        self.peak_current = 0
        self.diode_step = 0

        self.r_inner = r_inner
        self.length = length
        self.layers = layers
        self.wire_d = wire_d
        self.turns = round(length * layers / wire_d)

        self.layer_turns = []
        for l in range(layers):
            self.layer_turns.append(round(length / wire_d))
            # print("Layer", l, self.layer_turns[l], "Turns")

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
        self.voltage_list = []
        self.energy_list = []
        self.temperature_list = []

        # print("Wire Length", self.wire_length, "m")
        # print("Wire Mass", self.wire_mass, "kg")
        print("Turns", self.turns)
        print("Coil Length", self.length)
        print("Wire Diameter", self.wire_d)
        print("Layers", self.layers)
        # print("Resistance", self.wire_resistance, "Ohm")
        # print("Inductance", self.inductance, "H")

    def __str__(self):
        return f"Coil, {self.turns}({self.layers})"


class Stage:
    def __init__(self, coil, position, trigger_position, capacitance, initial_voltage):
        self.coil = coil
        self.pos = position  # Refers to axial distance of coil center
        self.trigger_pos = trigger_position

        self.capacitance = capacitance
        self.initial_voltage = initial_voltage
        self.voltage = 0

        coil.pos = position
        # coil.resistance = 2 * math.sqrt(coil.inductance / self.capacitance)
        # coil.added_resistance = coil.resistance - coil.wire_resistance
        # coil.tau = coil.inductance / coil.resistance
        # print("Added resistance", coil.added_resistance, "Ohms")

    def __str__(self):
        return f"Driver, {self}({self.pos})"


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

    def __str__(self):
        return f"Projectile, {self}({self.pos})"


class Material:
    def __init__(self, density, heat_capacity, resistivity, alpha):
        self.density = density
        self.heat_capacity = heat_capacity
        self.resistivity = resistivity
        self.alpha = alpha

