import math
import coil

Coil = coil.Coil
class Stage:
    def __init__(self, coil, position, trigger_position, capacitance, voltage):
        self.coil = coil
        self.pos = position  # Refers to axial distance of coil center
        self.trigger_pos = trigger_position

        self.capacitance = capacitance
        self.voltage = voltage

        coil.resistance = 2 * math.sqrt(coil.inductance / coil.capacitance)
        coil.added_resistance = coil.resistance - coil.wire_resistance


class Projectile:
    def __init__(self, coil, position, mass):
        self.coil = coil
        self.pos = position  # Refers to axial distance of coil center
        self.mass = mass

