// -*- mode: c++; -*-

/////////////////////////////////////////////////////////////////////

// IUPAC chemical element symbols

// References:
//
// J. Corish & M. Rosenblatt (IUPAC), Pure Appl. Chem. 75, 1613
// (2003).
//
// J. Corish & M. Rosenblatt (IUPAC), Pure Appl. Chem. 76, 2101
// (2004).
//
// W. H. Koppenol (IUPAC), Pure Appl. Chem. 74, 787 (2002).
//
// R. D. Loss (IUPAC), Pure Appl. Chem. 75, 1107 (2001).

// Number of official IUPAC elements
const int __niupac = 111;
// IUPAC element symbols
static const char *__iupac_symbol[__niupac] = {
	"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
	"Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
	"Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
	"Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
	"Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
	"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
	"Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
	"Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
	"Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
	"Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"
};
// IUPAC systematic symbol Latin abbreviation letters
static const char __iupac_systematic_symbol_rule[10] = {
	'n', 'u', 'b', 't', 'q', 'p', 'h', 's', 'o', 'e'
};
