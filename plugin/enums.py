import enum


class OverlayMethodEnum(enum.Enum):
    ALPHA_CARBONS_ONLY = enum.auto()
    HEAVY_ATOMS_ONLY = enum.auto()


class AlignmentModeEnum(enum.Enum):
    ENTRY = 'entry'
    CHAIN = 'chain'
    BINDING_SITE = 'binding_site'
    SELECTION = 'selection'
