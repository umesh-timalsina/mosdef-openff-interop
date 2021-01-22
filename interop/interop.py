from mbuild.lib.molecules import Ethane
from gmso import Topology as GMSOTopology
from gmso.external import from_mbuild
from openff.toolkit.topology import Topology as OpenFFTopology
from openff.toolkit.topology import Molecule, Atom


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
            openff_topology._box_vectors = gmso_toplogy.box.get_vectors()

        return openff_topology


def test_conversion():
    ethane_gmso = from_mbuild(Ethane())
    ethane_gmso.name = 'ethane'
    InterPackageConverter.convert_to_openff(ethane_gmso)


if __name__ == '__main__':
    test_conversion()
