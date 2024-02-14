import math
import numpy as np
import inputs


class Coil:
    def __init__(self, material, position, r_inner, length, layers, wire_d):
        self.type = type
        self.pos = position  # Refers to axial distance of coil center
        self.current = 0
        self.energy = 0
        self.temp = inputs.ambient_temp

        self.mass = 0
        self.r_inner = r_inner
        self.length = length

        self.layers = layers
        self.wire_d = wire_d
        self.turns = round(length * layers / wire_d)
        self.r_outer = r_inner + (layers - 1) * wire_d
        self.r_avg = (r_inner + self.r_outer) / 2
        self.area = pow(self.r_avg, 2) * np.pi

        self.wire_length = self.turns * np.pi * self.r_avg * 2
        self.wire_area = math.pi * pow(self.wire_d / 2, 2)
        self.wire_mass = self.wire_area * self.wire_length * material.density
        self.thermal_mass = self.wire_mass * material.heat_capacity
        self.wire_resistance = material.resistivity * self.wire_length / self.wire_area

        self.inductance = inputs.mu0 * self.area * pow(self.turns, 2) / length
        self.resistance = self.wire_resistance

        self.tau = self.inductance / self.resistance

        print("Wire Length", self.wire_length, "m")
        print("Turns", self.turns)
        print("Inductance", self.inductance, "H")
