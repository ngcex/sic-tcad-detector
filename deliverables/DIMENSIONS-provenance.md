# Layout dimensions: what is documented, what is not (2026-10-17)

Result of a read-only source review (papers/, DESIGN docs, public literature). No GDSII is
produced from invented numbers. The guard-ring numbers in the model (50/5/3/1 um, N_A 5e18,
etna/core/alternative_structures.py create_guard_ring_device, notebooks 19/20) are MESH
DEFAULTS, not device dimensions.

| item | value | class | source |
|---|---|---|---|
| active p+ area, dosimeter | 2 x 2 mm (4 mm2) | DOCUMENTED | JINST 20 C08019 p.3; Photons II.A p.3 |
| active p+ area, microdosimeter | 5 x 5 mm (25 mm2) | DOCUMENTED | Phys. Med. 140 (2025) 105675 s2.1; Altana, Sensors 23:6522 |
| epi / p+ | 10 um, N_D 0.5-1e14; p+ 0.3 um (0.3-1) ~1e19 | DOCUMENTED | JINST Fig.1; Phys. Med. s2.1 |
| substrate | 350 um n+ (Micro) / 100 um 'passive' (JINST Fig.1) / 'p+' (Photons text) | DOCUMENTED, CONFLICTING | three papers |
| front metal | Ni2Si ohmic + Al + Ti; Ti/Al bond metal periphery only | DOCUMENTED | Phys. Med. Fig.1a; Altana s3 |
| passivation | SiO2 + polyimide (7 um, SiCILIA 5x5) | DOCUMENTED (family) | Altana s5 |
| die size, 2 mm device | 2.7 +- 0.15 mm square | MEASURED FROM FIGURE | Photons Fig.1b (81/79 px = 2 mm, ~25 um/px) |
| periphery active edge -> die edge | 0.2-0.4 mm per side (+-0.1) | MEASURED FROM FIGURE | Photons Fig.1b |
| ring/frame lines | ~0.1 and ~0.2 mm inside die edge (+-0.05, resolution limit) | MEASURED FROM FIGURE | Photons Fig.1b |
| guard ring count | 1 drawn (JINST schematic, not to scale); none within 33 um of anode (Micro TCAD); JTE / floating rings / field plate are process options | QUALITATIVE, CONFLICTING | JINST Fig.1; Phys. Med. Fig.1b; Tudisco 2018 |
| field-plate overlap | 1-3 um on a 4.6 um dielectric (TCAD figure, +-0.3 um) | MEASURED FROM FIGURE | Phys. Med. Fig.1b |
| lateral depletion beyond anode edge | 10-15 um at 10 V | MEASURED FROM FIGURE | Phys. Med. Fig.1b |
| bond pad | 1 wire bond at a die corner; size unknown | QUALITATIVE | Photons Fig.1b |
| dicing street | unknown | MISSING | - |
| ring widths/gaps, oxide thickness, pad size, 5x5 die size | unknown | MISSING: owner/fab (SiCILIA mask set: STM / IMM-CNR via G. Petringa) | - |
| package | Al holder ~10 mm dia x 21 mm, 1 mm epoxy window | DOCUMENTED | Photons Fig.1d |

Other labs' typical designs (CNM, IHEP, LGAD) are not this device and are not used.

Recommendation: a GDSII may be produced only as a labelled schematic: top cell name containing
`REFERENCE-NOT-FAB`, every polygon labelled with class + source, one hatched 'termination zone'
0.2-0.4 mm instead of discrete rings, a placeholder corner pad marked SIZE UNKNOWN, no dicing
layer; DESIGN-2's 100 um SV array has no documented counterpart (the real microdosimeter is one
5x5 mm diode). Real dimensions require the mask data from the device owner.
