from mbuild.lib.molecules import Ethane
from openforcefield.topology import FrozenMolecule
from collections import OrderedDict
from gmso import Topology as GMSOTopology
from gmso import Atom
from gmso.external import from_mbuild
from openff.toolkit.topology import Topology as OpenFFTopology
from openff.toolkit.topology import Molecule
import simtk.unit as simtk_unit


class InterPackageConverter:
    @staticmethod
    def convert_to_openff(gmso_toplogy: GMSOTopology) -> OpenFFTopology:
        """Tries to convert a GMSO Topology into an OpenFF Topology"""
        molecule = Molecule()
        molecule.name = gmso_toplogy.name
        atom_indices = {}
        for atom in gmso_toplogy.sites:
            if atom.charge is None:
                charge = 0
            else:
                charge = atom.charge.value
            atom_indices[id(atom)] = molecule.add_atom(
                atomic_number=atom.element.atomic_number,
                formal_charge=charge,
                is_aromatic=False
            )

        for bond in gmso_toplogy.bonds:
            atom1 = atom_indices[id(bond.connection_members[0])]
            atom2 = atom_indices[id(bond.connection_members[1])]
            molecule.add_bond(
                atom1,
                atom2,
                bond_order=1,
                is_aromatic=False
            )

        openff_topology = OpenFFTopology.from_molecules(
            [molecule]
        )

        if gmso_toplogy.box is not None:
            openff_topology.box_vectors = gmso_toplogy.box.get_vectors().in_units(
                'angstrom').value * simtk_unit.angstroms

        return openff_topology

    @staticmethod
    def to_openff_frozenmolecule(topology=None):
        """
        Return an openforcefield toolkit FrozenMolecule object.

        This will return an OpenForceField Toolkit FrozenMolecule from a GMSO
        topology.

        Parameters
        ----------
        topology : `gmso.Topology` object
            An untyped Topology object

        Returns
        -------
        molecule : An `openff.toolkit.topology.FrozenMolecule` object
            A FrozenMolecule from the openforcefield toolkit

        NOTE:
        -----
        The current implementation converts a `Topology` to a `FrozenMolecule`
        using the `from_dict` method.
        """

        if not topology:
            raise TypeError(f'Expected: {type(GMSOTopology)}, '
                            f'Provided: {type(topology)}')
        mol_dict = OrderedDict()
        mol_dict['name'] = topology.name
        mol_dict['atoms'] = _create_atom_container(topology)
        mol_dict['bonds'] = _create_bond_container(topology)
        mol_dict['virtual_sites'] = tuple()
        mol_dict['partial_charges'] = None
        mol_dict['conformers'] = None
        mol_dict['properties'] = None

        return FrozenMolecule.from_dict(mol_dict)


def _create_atom_container(topology):
    if not isinstance(topology, GMSOTopology):
        raise TypeError(f'Expected {type(GMSOTopology)}, '
                        f'Provided: {type(topology)}')
    atom_tuple = tuple()
    for atom in topology.sites:
        if not isinstance(atom, Atom):
            raise TypeError(f'Expected {type(Atom)}, '
                            f'Provided: {type(atom)}')
        # NOTE: This is not the correct method to calculate "formal charge"
        charge = 0.0 if (atom.charge is None) else atom.charge.value
        charge *= simtk_unit.elementary_charge
        tmp_dict = {'atomic_number': atom.element.atomic_number,
                    'formal_charge': charge,
                    'is_aromatic': False}
        atom_tuple += (tmp_dict,)

    return atom_tuple


def _create_bond_container(topology):
    if not isinstance(topology, GMSOTopology):
        raise TypeError(f'Expected {type(GMSOTopology)}, '
                        f'Provided: {type(topology)}')

    bond_tuple = tuple()
    for bond in topology.bonds:
        tmp_dict = dict()
        tmp_dict['atom1'] = topology.get_index(bond.connection_members[0])
        tmp_dict['atom2'] = topology.get_index(bond.connection_members[1])
        tmp_dict['bond_order'] = 1
        tmp_dict['is_aromatic'] = False

        bond_tuple += (tmp_dict,)

    return bond_tuple


def test_conversion():
    ethane_gmso = from_mbuild(Ethane())
    ethane_gmso.name = 'ethane'
    InterPackageConverter.convert_to_openff(ethane_gmso)

def test_openff_conversion():
    ethane_gmso = from_mbuild(Ethane())
    ethane_gmso.name = 'ethane'
    InterPackageConverter.to_openff_frozenmolecule(ethane_gmso)


if __name__ == '__main__':
    test_conversion()
    test_openff_conversion()
