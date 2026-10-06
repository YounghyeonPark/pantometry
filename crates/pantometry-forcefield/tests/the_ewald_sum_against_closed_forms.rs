//! **The Ewald sum against what it must equal**: two Madelung constants, a lone charge's Wigner
//! lattice, the independence of α, the non-periodic limit of a neutral cluster at its stated
//! rate, the gradient, the strain derivative, and its own erfc against mpmath.
//!
//! Point charges here, with exclusions given by hand: no molecule and no force field, so that
//! what is checked is the sum alone. `tests/a_molecule_in_a_periodic_box.rs` checks it inside
//! UFF.

// The sums over atoms and axes are written as the formulas' index loops, which read against
// them; iterator chains over two arrays at once would not.
#![allow(clippy::needless_range_loop)]

use pantometry_forcefield::ewald::{erf, erfc, COULOMB, ERFC_ULP_BOUND};
use pantometry_forcefield::{Ewald, EwaldParameters, PeriodicBox};

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;

/// `(x, hi, lo)` bit patterns, `erfc x = hi + lo` from mpmath at 50 digits:
/// `python tools/erfc-reference/generate.py --points`.
#[rustfmt::skip]
const REFERENCE: [(u64, u64, u64); 438] = [
    (0x3fdc000000000000, 0x3fe127bf18c8eadc, 0x3c77ff0a3296d9cc),
    (0x3fdc000000000001, 0x3fe127bf18c8eadc, 0xbc81d1e2ab843044),
    (0x3fdbffffffffffff, 0x3fe127bf18c8eadd, 0xbc862f1321e4f5f0),
    (0xbfdc000000000000, 0x3ff76c20739b8a92, 0xbc77ff0a3296d9cc),
    (0xbfdc000000000001, 0x3ff76c20739b8a92, 0x3c81d1e2ab843044),
    (0xbfdbffffffffffff, 0x3ff76c20739b8a92, 0xbc94e8766f0d8508),
    (0x3ff0000000000000, 0x3fc4226162fbddd5, 0xbc4b40443f6ec34a),
    (0x3ff0000000000001, 0x3fc4226162fbddd2, 0xbc6b59095b5991f3),
    (0x3fefffffffffffff, 0x3fc4226162fbddd7, 0xbc6c8b94ea1cc038),
    (0xbff0000000000000, 0x3ffd7bb3d3a08445, 0x3c98da0221fb761a),
    (0xbff0000000000001, 0x3ffd7bb3d3a08446, 0xbc8929bda9299b83),
    (0xbfefffffffffffff, 0x3ffd7bb3d3a08445, 0x3c8722e53a87300e),
    (0x3ff8000000000000, 0x3fa15aaa8ec85205, 0xbc2e86ee834da4ce),
    (0x3ff8000000000001, 0x3fa15aaa8ec85201, 0x3c23322ca5adc672),
    (0x3ff7ffffffffffff, 0x3fa15aaa8ec85209, 0xbc4410026b1243d6),
    (0xbff8000000000000, 0x3fff752aab89bd70, 0xbc8385e445f2c96d),
    (0xbff8000000000001, 0x3fff752aab89bd70, 0xbc613322ca5adc67),
    (0xbff7ffffffffffff, 0x3fff752aab89bd70, 0xbc915f7feca76de1),
    (0x4000000000000000, 0x3f7328f5ec350e67, 0xbc1ca006412e68d0),
    (0x4000000000000001, 0x3f7328f5ec350e5c, 0xbbdd7435b21abed4),
    (0x3fffffffffffffff, 0x3f7328f5ec350e6c, 0xbc0408cf68698d84),
    (0xc000000000000000, 0x3fffecd70a13caf2, 0xbc99a35ff9bed197),
    (0xc000000000000001, 0x3fffecd70a13caf2, 0xbc96fe28bca4de54),
    (0xbfffffffffffffff, 0x3fffecd70a13caf2, 0xbc9af5fb984bcb39),
    (0x4004000000000000, 0x3f3aab859b20ac9e, 0x3bd88f4ff748376b),
    (0x4004000000000001, 0x3f3aab859b20ac8d, 0xbbdd7d0c47b5860f),
    (0x4003ffffffffffff, 0x3f3aab859b20acb0, 0x3bcd37586c8bef5d),
    (0xc004000000000000, 0x3ffffe5547a64df5, 0x3c8b0cee160116f9),
    (0xc004000000000001, 0x3ffffe5547a64df5, 0x3c8b9bafa188f6b1),
    (0xc003ffffffffffff, 0x3ffffe5547a64df5, 0x3c8a7e2c8a793741),
    (0x4008000000000000, 0x3ef729df6503422a, 0x3b6784ca4c429a15),
    (0x4008000000000001, 0x3ef729df65034218, 0xbb8a661b3e88955f),
    (0x4007ffffffffffff, 0x3ef729df6503423c, 0x3b9314403254f4a1),
    (0xc008000000000000, 0x3fffffe8d6209afd, 0xbc908a82f0994988),
    (0xc008000000000001, 0x3fffffe8d6209afd, 0xbc9085f2ccf260bc),
    (0xc007ffffffffffff, 0x3fffffe8d6209afd, 0xbc908f1314403255),
    (0x400c000000000000, 0x3ea8ef2a9a18d857, 0xbb42d76dc03e80a5),
    (0x400c000000000001, 0x3ea8ef2a9a18d840, 0x3b1e164ec8e45b8c),
    (0x400bffffffffffff, 0x3ea8ef2a9a18d86d, 0x3b468e5aa6667838),
    (0xc00c000000000000, 0x3fffffff3886ab2f, 0x3c8c9ea52d76dc04),
    (0xc00c000000000001, 0x3fffffff3886ab2f, 0x3c8c9effc3d3626e),
    (0xc00bffffffffffff, 0x3fffffff3886ab2f, 0x3c8c9e4a971a559a),
    (0x4010000000000000, 0x3e508ddd13bd35e7, 0xbaf615db40319381),
    (0x4010000000000001, 0x3e508ddd13bd35c5, 0xbafba01f978aa578),
    (0x400fffffffffffff, 0x3e508ddd13bd35f8, 0xbaf350b914850422),
    (0xc010000000000000, 0x3ffffffffbdc88bb, 0x3c70b2865615db40),
    (0xc010000000000001, 0x3ffffffffbdc88bb, 0x3c70b28edba01f98),
    (0xc00fffffffffffff, 0x3ffffffffbdc88bb, 0x3c70b2821350b915),
    (0x4014000000000000, 0x3d7b0c1a759f7739, 0xba1b28568855c7a1),
    (0x4014000000000001, 0x3d7b0c1a759f76f4, 0xba161d46d7ef49d6),
    (0x4013ffffffffffff, 0x3d7b0c1a759f777d, 0x3a1fcc99c743e5a7),
    (0xc014000000000000, 0x3fffffffffffe4f4, 0xbc7a759f7738935f),
    (0xc014000000000001, 0x3fffffffffffe4f4, 0xbc7a759f76f3a78b),
    (0xc013ffffffffffff, 0x3fffffffffffe4f4, 0xbc7a759f777d7f32),
    (0x4018000000000000, 0x3c78cf81557d20b6, 0x38fa7fff0cc732c0),
    (0x4018000000000001, 0x3c78cf81557d206b, 0xb9157190d240b5c1),
    (0x4017ffffffffffff, 0x3c78cf81557d2102, 0xb91d4e6fa75b784b),
    (0xc018000000000000, 0x4000000000000000, 0xbc78cf81557d20b6),
    (0xc018000000000001, 0x4000000000000000, 0xbc78cf81557d206b),
    (0xc017ffffffffffff, 0x4000000000000000, 0xbc78cf81557d2102),
    (0x3fdc000000000000, 0x3fe127bf18c8eadc, 0x3c77ff0a3296d9cc),
    (0x3fdc000000000001, 0x3fe127bf18c8eadc, 0xbc81d1e2ab843044),
    (0x3fdbffffffffffff, 0x3fe127bf18c8eadd, 0xbc862f1321e4f5f0),
    (0xbfdc000000000000, 0x3ff76c20739b8a92, 0xbc77ff0a3296d9cc),
    (0xbfdc000000000001, 0x3ff76c20739b8a92, 0x3c81d1e2ab843044),
    (0xbfdbffffffffffff, 0x3ff76c20739b8a92, 0xbc94e8766f0d8508),
    (0x4018000000000000, 0x3c78cf81557d20b6, 0x38fa7fff0cc732c0),
    (0x4018000000000001, 0x3c78cf81557d206b, 0xb9157190d240b5c1),
    (0x4017ffffffffffff, 0x3c78cf81557d2102, 0xb91d4e6fa75b784b),
    (0xc018000000000000, 0x4000000000000000, 0xbc78cf81557d20b6),
    (0xc018000000000001, 0x4000000000000000, 0xbc78cf81557d206b),
    (0xc017ffffffffffff, 0x4000000000000000, 0xbc78cf81557d2102),
    (0x0000000000000000, 0x3ff0000000000000, 0x0000000000000000),
    (0x0000000000000001, 0x3ff0000000000000, 0x0000000000000000),
    (0x8000000000000000, 0x3ff0000000000000, 0x0000000000000000),
    (0x8000000000000001, 0x3ff0000000000000, 0x0000000000000000),
    (0x403a800000000000, 0x0043df6725a60cf5, 0x0000000000000001),
    (0x403a800000000001, 0x0043df6725a608d7, 0x0000000000000001),
    (0x403a7fffffffffff, 0x0043df6725a61113, 0x0000000000000001),
    (0xc03a800000000000, 0x4000000000000000, 0x0000000000000000),
    (0xc03a800000000001, 0x4000000000000000, 0x0000000000000000),
    (0xc03a7fffffffffff, 0x4000000000000000, 0x0000000000000000),
    (0xc03b000000000000, 0x4000000000000000, 0x0000000000000000),
    (0xc03b000000000001, 0x4000000000000000, 0x0000000000000000),
    (0xc03affffffffffff, 0x4000000000000000, 0x0000000000000000),
    (0x403390868aaf9bd1, 0x1d1944fc6842c4af, 0x99772e2541d8520c),
    (0x40399d96f5d4243b, 0x046d005cfa377d18, 0x00c615b498e35202),
    (0x40367bff222b2f1b, 0x120425b57aa6b814, 0x0e9b778f14dfe10a),
    (0x40353965e0bc0ccd, 0x16fd69764511f626, 0x139a20d4e06d76bf),
    (0x4037296acb8ce60b, 0x0f399e0909cd900f, 0x0bd7c21ece580d74),
    (0x40375b392ec5a402, 0x0e686bfbb8703683, 0x0af92fba1602c320),
    (0x40375bfe647c85a2, 0x0e6537626f1b7ed0, 0x8b00cd1866d3a34f),
    (0x402f1088e289f5dc, 0x29e1e2dbb8900487, 0xa679019c64ea8c67),
    (0xc0007e125e7756a5, 0x3ffff1741d61f0c8, 0xbc8f1cdcef571262),
    (0x403067a4ce076348, 0x275d697c4879dfb2, 0xa3f55f5191761779),
    (0x400ea2a449db7254, 0x3e706561855135fa, 0xbadfe4d0070477c1),
    (0x402d08f01aa74120, 0x2ca3137e91f22728, 0x29469a7d879d7263),
    (0x40398ab96e783862, 0x04c3b6818dbc6e58, 0x81605310defd3026),
    (0x40291c680ed155df, 0x3171104d4395ccdc, 0x2e07c2818ae519b9),
    (0x402d9b2216891418, 0x2be1b876e6c5a41e, 0x288904a83ee08e1f),
    (0x40375136ff0a715f, 0x0e92f6578b6c7581, 0x8b076d8449406b2e),
    (0xbff29d16f5c58b8c, 0x3ffe66b761023721, 0x3c97bfd9d91d0135),
    (0x402039ea10db8f62, 0x39c22464b56b51ea, 0x366ea7219827bbbe),
    (0x40157657575c6268, 0x3d2247eaffa5872f, 0x39bd3ce1f23d88e4),
    (0x4033c722807a2467, 0x1c57330fc08891a1, 0x18f79d915b366e28),
    (0x4025cb931c4ffed6, 0x34f4faf647b375ed, 0x318c97aea493c2d9),
    (0x401755a938bb1da0, 0x3ca6d22d8d05ad34, 0x394807082979d988),
    (0x4030138b4d5eb2be, 0x2853b8d7725e0b60, 0x24ecf580a837711d),
    (0x4009d2e985e974bc, 0x3ed4f1c2bec95b92, 0xbb7417ff17c77664),
    (0x402dc286f0176580, 0x2bacbd9740abc15e, 0xa83e5ad0204d8047),
    (0x3ff46a4044b31630, 0x3fb23792e16ae126, 0x3c550ffb8caa1ef6),
    (0x400ec8e66f0b1a8c, 0x3e6c47e0422ffc08, 0x3af9a3b01e8eab7c),
    (0x401e9ddde4eee248, 0x3a6a0fa4db2ecfd4, 0xb70af7d813a05076),
    (0x403857556abfed28, 0x0a2b8935917994f4, 0x06ca90c4e08afb36),
    (0x401fd8c74d8cc3da, 0x39fa58be47cf9386, 0xb69c0f57c0bf92fc),
    (0x4018fa28aaf1b16a, 0x3c332d5b0c86e214, 0x38d2aaf0f1c4d8f2),
    (0x4025315e8f5903bc, 0x358b4b4635c9c9e9, 0x322cf683de2729a0),
    (0x4003abdca5f11dda, 0x3f4096a258ae4614, 0x3bb2e38703d99c5c),
    (0xc001826ed7c02442, 0x3ffff7f229172cbb, 0xbc65ab0a321c45cd),
    (0xbfd3b6c468451c68, 0x3ff563e496e46789, 0x3c8856e3f692a0f7),
    (0x403a16711d60eb43, 0x023890380b64c822, 0x000000013cebada4),
    (0x4023a774b51677d2, 0x36f762113cd60ad0, 0x339e469171919ec3),
    (0x402fc2166344916c, 0x28e550f5c56e4ca0, 0xa56cffcbf39750a5),
    (0x402c8234c3eece92, 0x2d526aadc762e4a3, 0xa9f5510ff9eb0375),
    (0x4033db9ec5a7f6d8, 0x1c0f012d5587b8bb, 0x9895269d6f4be91b),
    (0xbff089704bf6f278, 0x3ffdb2dc9787d472, 0xbc7ba7e1c86c0554),
    (0x40337c3b9fa425fa, 0x1d6185037aee95b4, 0x1a09badf913588b1),
    (0xbfe390dea8b29cf4, 0x3ff9ce02f8e67227, 0x3c9a87148dea4e1e),
    (0x4030bfb3c6f92ead, 0x26550cefb8c287e5, 0x22fad646886aaf16),
    (0x40360c16430a93d5, 0x13c5c9e9d39b0734, 0x1056d96200718db1),
    (0x4032e0a488b6673e, 0x1f7c2253e78287ef, 0x9bf46fb87c4a0a06),
    (0x4014b05a4cbc1920, 0x3d52299bee0832cb, 0xb9e0c91cab19c1e4),
    (0x4024af5d3de2ec06, 0x36063f8ffe32874b, 0xb252673e91d081d3),
    (0x402dc46b1194fa94, 0x2ba9bdfa5cf07e49, 0x284936c8908938ef),
    (0xbfd636170a0166a0, 0x3ff605e134315c96, 0x3c92637f74aafa69),
    (0x4026aa3c94932b56, 0x3414ee6ff832645c, 0xb0b0e9efd1ada3b8),
    (0x400a72d52fac000c, 0x3ec8998d0e646347, 0x3b56cad52cf5f138),
    (0x4015024546095396, 0x3d3f1717808c5133, 0x399a9cd4d4b60d6a),
    (0x4029ab544d064c4a, 0x30cc9083f5375c00, 0xad5d1da0642d2ddb),
    (0x3fee521851a50d20, 0x3fc7124e320f1bf7, 0x3c552deb3c015a54),
    (0x4004269416f4f5d4, 0x3f381a9bed15393a, 0xbbdb3c6f0a60dc5c),
    (0x40026b29d537aacc, 0x3f5283a08a932dc8, 0x3bd6f2daff2b6625),
    (0x4012189b2c831736, 0x3de5a4f2dc58a7a4, 0x3a49b0583f67e02e),
    (0x401537865cff4e46, 0x3d31c901e6c80342, 0xb9d43c0fbd541cf2),
    (0x3ffe2e6640de1338, 0x3f7f48edc3fe58f4, 0xbc1d5afbb800075a),
    (0x403a593e3014956a, 0x00fd1eee82861fb9, 0x0000000000001400),
    (0x40337938cd12ab7c, 0x1d6bb846b92dcd00, 0x19ea938de15f90d7),
    (0x4036f8c56468a619, 0x10043a49dee0b14e, 0x8c9839fe36b85942),
    (0x4033e82150be430c, 0x1be1b8cf270073ee, 0x1879f44057da8f4d),
    (0x402bb383e2073800, 0x2e58772173c11a30, 0xaae871a3fa40ba2c),
    (0x402cbfe58e9d0704, 0x2d028fdab451f659, 0xa976ce37c5c2ecb3),
    (0x3fc820848cf921b0, 0x3fe9460dadd1667b, 0xbc7db006dd989e0c),
    (0x4023110b45441f0d, 0x377bc86d58232a59, 0x341637ef84f2232f),
    (0x3feaeca13c1f4844, 0x3fcdf68947293ec8, 0xbc6d04a31d8a78c6),
    (0xbfc0f8fd3c598640, 0x3ff26146c87150cf, 0x3c9393b2e7c6a75a),
    (0x402f07c8ba1681f6, 0x29ee70a54fe789e4, 0x266ee1ba0a45c36c),
    (0x40281b85b3422c64, 0x328f3884144f5b97, 0x2f174af8d3a34018),
    (0xc0009ac86b8540fa, 0x3ffff259fdc1ae19, 0xbc795c29cd3ab19b),
    (0x4033dc747f8c4c2a, 0x1c0b3bb879221c2c, 0x98ae3357fc38cdea),
    (0x40364f7e8c462084, 0x12b7e24dd1c8e7d4, 0x0f5af4a1d43bab20),
    (0x401ccbf0e46c8d1a, 0x3b07475bd7cb3857, 0xb76c1c7c50fca64b),
    (0x40313d8082ba5011, 0x24d32123278b9a41, 0xa17fe6ebf43d9b81),
    (0x401761cec2053e7c, 0x3ca3d456daccb6fd, 0xb8fcb254bf538baf),
    (0x40285a1a330aa1ec, 0x32499691d6f18b9a, 0xaedc7b50c055621d),
    (0x3fbae30059723200, 0x3fec38ba77588bd4, 0x3c8cd2a506b7ead6),
    (0x401dcd5cab77d816, 0x3ab21dc1b5e4114a, 0x375a9e06f2834c87),
    (0x4030c11b7f3a9c32, 0x265182417824b2bf, 0x22fd389b85f28daf),
    (0x3f93f05b07a71e00, 0x3fef4c08a27348cf, 0x3c894a68e6e1e3a5),
    (0x403382774ae539b6, 0x1d4b155983c2ee08, 0x99bd2136fed6fd6a),
    (0x4021235552c6e670, 0x3911866397072dd6, 0xb5bbf1e80c9e2c98),
    (0x4032aa9a9e530386, 0x20333642de664c1f, 0x1ca8fb22fb139bca),
    (0xbfe0d9583e545fe0, 0x3ff8b2306da65f9a, 0x3c3cacb2be9dccc9),
    (0x402abd1f530115d6, 0x2f879d674c5339be, 0x2c2a808cccaf4921),
    (0x4001ffaa7718ceec, 0x3f57fbf797ad5152, 0x3bbdca062cb18847),
    (0x402f447ed28ed938, 0x29980e383eb22fe9, 0xa636528afc2efad1),
    (0x40363997a0db6691, 0x1310eb0a487b8980, 0x8fb3782d1d73ba22),
    (0x4036d31338936f75, 0x10a0de8c55977e10, 0x8d4d20c00b1d6c55),
    (0xbfdfb4bdae250d70, 0x3ff843653f351521, 0x3c4a730a2ed99d10),
    (0x40117ab372f8da1c, 0x3e060b39423940e3, 0x3a92918dce71285d),
    (0x4036f15dfd2630c7, 0x10231ad17029a91f, 0x0cc7dc69e2a53479),
    (0x403497f8e9d88963, 0x195f375ab3b593f4, 0x95f99487b53bc042),
    (0x403132cc9778181a, 0x24f43ad9ca571fa3, 0x2188789d250edbff),
    (0x400b8321b374d900, 0x3eb358c27059e3ce, 0xbb50ee173232480e),
    (0x402e9ad839b0fb9e, 0x2a8539989cfe7848, 0x27257df090461d25),
    (0x400463f0d4974816, 0x3f347e36442a64d7, 0x3bdc246954ce5960),
    (0x4035affcca40f254, 0x1531fc4dfe23cc37, 0x11d97df8d3fc58f9),
    (0x3ff0bd251be3adc8, 0x3fc1cace41065dc9, 0xbc6cab58ffe388b8),
    (0x402483317f82a189, 0x363094804cb338f2, 0x32d11461ea5346ce),
    (0x40324bc11638e042, 0x217090fe6ca673af, 0x1e185a5e1c38b1b3),
    (0x402fbf71337dba86, 0x28e91fe415bbff56, 0x257ffaf887551da4),
    (0xbff40ccf1080851b, 0x3ffec7377773b0e1, 0xbc9e37b02cb06cd4),
    (0x402506bc10a7e6d4, 0x35b3f153cccb9903, 0x323dc2ebb0c63fc3),
    (0x4020a843c1bf928f, 0x39705a9e04a4fe48, 0x361f89a8e3531363),
    (0x4036a5d05333a277, 0x1159be94beb36413, 0x0dfd3f3f4ed888d8),
    (0x40369fdfec65324f, 0x11726a9783ad2252, 0x0e1719a128d11e13),
    (0x4033548a932514be, 0x1dec570a70a55b9f, 0x9a53fb622b866f52),
    (0xbff62bccf132ac60, 0x3fff330ebb3bb973, 0x3c8c73a5196c24bb),
    (0x40183c47cea4bb9e, 0x3c682aff3fae65bf, 0xb8ea1e5b988fbcd1),
    (0x4029df50a840f5cd, 0x309090034f7f7905, 0x2d176ce297d68a65),
    (0x401ee1bdc7a010de, 0x3a52a73cfae437d2, 0x36ffadae378f71b4),
    (0x4034cd06973257a9, 0x1897cf6bac588121, 0x153d5276a92f7331),
    (0x403636c461ca27fc, 0x131ba41d8f483076, 0x0faac2c6d0fe6b06),
    (0x40386516b5cb210a, 0x09f00640ed00715e, 0x0688be91973cd01a),
    (0x402f9912730901ec, 0x2920f0416f7d07ce, 0x25c8d64d32e2a92d),
    (0x4031c992ef4b9ca1, 0x2317a1e500c31801, 0x1fb409aeab41f990),
    (0x3fe2484e229e1a78, 0x3fdad297ccc04a13, 0xbc7abdc04559e730),
    (0x4027706525dac018, 0x334621edf4c7fbec, 0xafeb8fe1d5a90987),
    (0x402e83b87cf4f82c, 0x2aa5287e4af47b24, 0x273aed3eb2989edc),
    (0x40171d0fa71709ca, 0x3cb5e1947535a11a, 0x395aef27d48c3188),
    (0x401260f5f42610ee, 0x3dd6649592c59e0e, 0xba7e866a70ce0ab8),
    (0x401e69cc5e287c8a, 0x3a7c7e44a1d4a35c, 0xb70763836a4f25af),
    (0x40398dd592149f27, 0x04b52effe5d03fa5, 0x015f8c5e0c38620b),
    (0x4026d8940873f91b, 0x33e52eb7c5e9aa31, 0x3072049731858ed5),
    (0x401708de343c36ec, 0x3cbb9003fa351958, 0x390be9c146986c7e),
    (0x4031b9950ac4e856, 0x234b4010932e7b06, 0x9fe7ad941336781e),
    (0x40032e02ec65f37c, 0x3f46dbb354e01aef, 0xbbe126c00c18569e),
    (0x40135a476ca9fb26, 0x3da1265ab99d9b79, 0xba22c0c673dd190e),
    (0x4037354dd839139a, 0x0f07c01d51f2ce9d, 0x0ba56340874e3e4e),
    (0x4033559e2e26c31b, 0x1de814f1597fc381, 0x1a6eb2aba31e3e8b),
    (0xc00669acdd5bc8ad, 0x3fffffb215a4a38e, 0x3c90967e89700718),
    (0x40302da0519b4d40, 0x280773dc658339b9, 0x24aa2d9bf623bedd),
    (0x40348a73f33c90d2, 0x19912bcc6f96f144, 0x962730a8485a82ae),
    (0x403037804984ba37, 0x27ead300011d7af3, 0x248b8621907d3e99),
    (0x40257c6a5461b12d, 0x3542dff4e3a43be9, 0xb1d9cd7885f07a8b),
    (0x403330f383a63573, 0x1e679c4baa78b282, 0x1afb556a1f6c3a47),
    (0x401e2b005b0eac3e, 0x3a922d2796c4d7b2, 0x373a21576412c1f1),
    (0xbfdf6231eed4f248, 0x3ff83121cc2dc021, 0xbc67217452a390ee),
    (0x40396c20964ce945, 0x055118acc554dcd7, 0x81d6406481f0da59),
    (0x4037d0d57e9f0ae5, 0x0c73ecb91c46ea81, 0x88e0454933bb46d1),
    (0x402f951e58db0678, 0x2925a0de7741a93e, 0x25bf0b2fc126caf6),
    (0x402790490f1e0f23, 0x33245fde71692162, 0x2faecf1ebca48437),
    (0x4031856447f60263, 0x23f1c541a6345f6a, 0x209696631c713b23),
    (0x402c823358cbd6f4, 0x2d526c23214c6f48, 0xa9eb2baa3a1c3c04),
    (0x402b19b3ea64966e, 0x2f16f0d5bdf281a8, 0xaba7553a1aa7b08b),
    (0x40001e6c2c9504ac, 0x3f71f071c7771d17, 0xbc1421f399ef0dcf),
    (0x40366062b4e593e5, 0x1273f8db9cd53a99, 0x8f1c86210ca0856f),
    (0x4037735d4d780892, 0x0e02d8fc1699adff, 0x0aa92471a1e2d834),
    (0x4037623bd0dd9f94, 0x0e4b21850f3b86d6, 0x0ae6a16ce05cc66b),
    (0x4024adef56dae92b, 0x360793b2e41bceb2, 0xb29f91048ea7a24e),
    (0x4038ed26459dcd33, 0x0791c9ed35a7ebf1, 0x04236fdbf6680b21),
    (0x4024bf2b18645956, 0x35f768104d5e8af6, 0x3292e2d0852db58b),
    (0x40124cca134f6c00, 0x3ddaefc753fb6c97, 0xba786a3a10ebc534),
    (0x401be244ff0edd14, 0x3b530c3deff63ef8, 0xb7eda36defbdcdc8),
    (0x401504960b2aa7f6, 0x3d3e58fe6ced3c35, 0xb9a4a8b5ba2488ad),
    (0x403709613259d19e, 0x0fc04ff4f349bc0e, 0x8c65f58a30008e00),
    (0x4034a78d5466c4d8, 0x19243abd7e16acdc, 0x95b092f9bd1ea1cd),
    (0xbff378b311d35ab1, 0x3ffea2d9f4931854, 0x3c83b7adc8faf592),
    (0x40209c89fe5189ad, 0x3978008876644e8b, 0xb60297ce9d6d040f),
    (0xbfc7dd0de03abe80, 0x3ff353c94c0bf566, 0x3c9857693103d8ce),
    (0x4017f8370913b90c, 0x3c7b36682258e9fe, 0x38f38397ac786324),
    (0x40308de16871d750, 0x26eb27b0cbeb21d0, 0x23818649933e02ee),
    (0x40372daa83d0afa2, 0x0f27b9e860e3d927, 0x8ba6fdc4e3b90295),
    (0xc006a74e36728a7e, 0x3fffffbedd58a97f, 0xbc93dc07417d181d),
    (0x401d8f6248f7272e, 0x3ac66b4385d58b36, 0xb73f512e96591032),
    (0x4016f0aa80848c44, 0x3cc227bc43726887, 0x39418bb3ddd01c4b),
    (0xc001a4a8530df19c, 0x3ffff89094a58774, 0xbc9c5851189a2737),
    (0xbfc7e8e802ef05b0, 0x3ff355669b5f71d7, 0xbc8db142eea6b9aa),
    (0x402d18ae6103f072, 0x2c8f2715cb98acf2, 0xa924f67f3b780ab7),
    (0x3ff1671c0b08cb4c, 0x3fbfbe68f4489610, 0x3c26f7c0eea5500c),
    (0x4035a657e721a3f1, 0x15570c87bfad1c0e, 0x11e8e26c207c7f55),
    (0x4035075d94b3fcff, 0x17bbf595f22de520, 0x944ee158164d8a21),
    (0x40314c56eb9e1522, 0x24a49abd0370569c, 0xa1434e03d592becc),
    (0x4026fb5a997d24ee, 0x33c1c3e6dd23f7bd, 0x306c99764471e06d),
    (0x4030ee35364e519e, 0x25c77b59a98c8c92, 0xa26735ded94ce4b9),
    (0x400492f70bec1ac6, 0x3f32135d3aec2d00, 0x3bc385c93ad278e8),
    (0x40390f3da6c35380, 0x06f74ae733928799, 0x839287790cdc63d2),
    (0x4034b94f85cd9346, 0x18e2488951790e0e, 0x155849212bb2b8ef),
    (0x4015234f2b6107be, 0x3d35ffd46e1bcdcf, 0x39d3fe3107510783),
    (0x4025441f663b5600, 0x3578ff97e27363b6, 0x31e5e269d387a83f),
    (0x4028f29d90ebe730, 0x31a08f3d9ec2c58d, 0x2e398dd819259de1),
    (0x400d0fa94039b87e, 0x3e92b34805ceccec, 0x3b30e4147d22ec75),
    (0x402c19b7aef9d5e0, 0x2dd78ad76a40bfdb, 0x2a2817b7bc6b14a9),
    (0x402ff9fa8d8d4ea2, 0x2894e5ed828234f5, 0xa53302feff9870c2),
    (0xc00605e824ed3ee1, 0x3fffff983fd512cc, 0xbc6e6a5f87342420),
    (0xc00125c7623fdbfe, 0x3ffff6069c8071fa, 0xbc9dc9808c71f8f2),
    (0x400992b522cdecde, 0x3ed9d94fac01a19d, 0x3b509dc38eb93360),
    (0x401b0c33eec95ed0, 0x3b95a6a453b8f7cf, 0xb83308fd49c875a2),
    (0x401f1c7a6629edb0, 0x3a3e72eacfb13a5f, 0xb6bf7ede18aff08a),
    (0x40386cee6238af3d, 0x09ccb1ccce4e1e17, 0x866394969ea3483a),
    (0x4033c15415fa7d83, 0x1c6c77beb8994da7, 0x1902bc881c21a463),
    (0x4039c53c6f5cb75a, 0x03b4a72585bc5e1a, 0x805b1156227681c6),
    (0x4030f4ddbd534434, 0x25b36c618d64988e, 0x225866d4e879235c),
    (0x401bd35b272b06ee, 0x3b57610900e172c2, 0xb7fb8e0cfe482a97),
    (0x401e24397b83e8b6, 0x3a9419d9d661e0b1, 0xb7333f94b99cebbe),
    (0xbffafe7ff58c49da, 0x3fffba3a48c4ce94, 0xbc83b7f46037a76b),
    (0x40245583db63121c, 0x3659dbb9ebd9966c, 0xb2f31a0a5aea2214),
    (0xbfd402a4a115da00, 0x3ff57754eb811310, 0xbc8275b297b1edec),
    (0x4016aacedeb258d0, 0x3cd3fda8baecbded, 0x396c3e6b4834561b),
    (0x4034d95304f56e10, 0x1869b24333230cf9, 0x94eeb0394f152272),
    (0x4038f6ce554328e4, 0x0765a4258abf0ac7, 0x840eb86e924f17fc),
    (0x402158a093e148d3, 0x38e7057946dce3d7, 0x3577967f026230a4),
    (0x403487a70004ab2a, 0x199aecdb9adad27c, 0x16270d2f35465668),
    (0x402e5869acb090c6, 0x2ae1714927f1c452, 0x278c3efe6d14f85b),
    (0x3fb617855b509c20, 0x3fece4477d0f2c58, 0x3c841c9d5808add2),
    (0x400295538413977a, 0x3f50b400da63512f, 0x3bd8a781256256c2),
    (0x4034c2abe28f4e44, 0x18c004427c2b5084, 0x9562faecc1da9b2d),
    (0x4015fc97b2627d38, 0x3d01322c5d13392c, 0x3935e3acbf1e0a5c),
    (0x402b04528c957f32, 0x2f31ceb6ccabdce0, 0x2bb3d4e7da7bceab),
    (0xc0038137d1c0ab95, 0x3ffffdafc8811f1c, 0x3c90110bc5ae3fa8),
    (0x40346828fc154bbb, 0x1a10470a0d369b4b, 0x96a3e457d14c3f54),
    (0x402bfbf117b6ce60, 0x2dfe3034f350776d, 0xaa71c1280fe9530e),
    (0x4036482c82d351f0, 0x12d5661e257d13c8, 0x0f581efcd3802149),
    (0x4021c3e5250be7f5, 0x38922f1615037f44, 0xb53b3c79d18ed04b),
    (0x4031dff7051f47ae, 0x22d0a0f18cdc446a, 0x9f7157bf2d3647b5),
    (0x3ff01eb3c0b50670, 0x3fc3bd2f6707df4a, 0xbc426f7b2cdf3c20),
    (0xbff8d73d71b89da4, 0x3fff8cd41a766dcd, 0x3c93cef6c7c251d8),
    (0x4039ee26c86609da, 0x02f5adc5f75a79e9, 0x8000153e509c3f83),
    (0x40388e01ff325bff, 0x093a15020429f6d9, 0x85946cc833d980c9),
    (0x400cc463aaa8f9e4, 0x3e989fc465b18f4b, 0xbb35c545dad62634),
    (0x3ffcab150a997850, 0x3f87194b859000b6, 0x3c07265dadb4ea6e),
    (0x4039a67c61790f3c, 0x044381b993b650cf, 0x00bb524a2544dbc5),
    (0x40207cfc62823bfb, 0x3990c30e912d0ff6, 0x363cd85737bad730),
    (0x4021054de43201e3, 0x3928059f4a7695f7, 0xb5c1f7eed62674fc),
    (0x4023e606d47a57be, 0x36c07bdd007f68a7, 0xb36332dfeded8597),
    (0x4009596815b9ee7c, 0x3edf2373d60f0324, 0x3b749d3aba144e09),
    (0x4036eb8155802c24, 0x103b544fdf88f6b7, 0x0cb920c83ca773e9),
    (0x40361b6de6e326db, 0x1388a8293412ba55, 0x100da6b3ce5c2792),
    (0x40200ba885e487cd, 0x39e3b3e168da2c76, 0x3679c19501093edb),
    (0x3ffd14103f32f223, 0x3f84d10489d93011, 0xbc2d5fd89859c79e),
    (0x400b54e2c992bdc2, 0x3eb6b97adc88a7dd, 0xbafca9dc2bc56326),
    (0x3fea88b1954eb184, 0x3fced70b45bf5e0c, 0xbc446ad07a19f141),
    (0x40016f059213fe57, 0x3f60d9f76e68214d, 0xbbfc687d5dd5790e),
    (0x401245771dd11a98, 0x3ddcccd15ea20067, 0x3a4a9fad1a0f8324),
    (0x3fd78666eed0ee84, 0x3fe34d3e6a13519f, 0xbc82ea8c79beeb2d),
    (0x3ff0e529ef2ca4e6, 0x3fc1531f4a3ef10d, 0xbc67bc84fcd38164),
    (0x3ff88719a6cb0ae8, 0x3f9ee298e481a41a, 0x3c15fd8099c4bce6),
    (0x400bf2141097f820, 0x3eaa321047e84f79, 0x3b42a6ec9e21abd0),
    (0x40006e34fea8a2aa, 0x3f6e217551dd1664, 0x3be00a36390808ce),
    (0x40003850db97e7b1, 0x3f70f4cbcf8e4853, 0xbc187d315c519664),
    (0x40115208b34f0f66, 0x3e0f69e47e49df95, 0xba8e7a31ca5fd932),
    (0x3fe2b352c237d3ab, 0x3fda25a49d1c99d9, 0xbc6ee0b8f2983e34),
    (0x3ffb2c73149ad07c, 0x3f90b46730f0413b, 0x3c3f81d4c2c48bef),
    (0x3fec98419b5571cf, 0x3fca6905d60d38dc, 0xbc620fd19e263060),
    (0x3fd50b19dedf56ac, 0x3fe48ab84b88a592, 0xbc42c85111280450),
    (0x400511fb6c9d7222, 0x3f29a0c3ca70cc3b, 0xbb6317c76c4b1410),
    (0x40029ef95c2c42e3, 0x3f504ffe1761b5e7, 0xbbe3d9fa0ceb1f15),
    (0x400b1641a92499bd, 0x3ebc36b6471e68bf, 0x3b5381a46d8efa5b),
    (0x4004815d76a7485f, 0x3f32f280f484e39f, 0xbb998335136ac210),
    (0x3ff3dc9d40c1529b, 0x3fb4442e932a9cbd, 0xbc5ef77d51a56c8a),
    (0x400233e6f62d9536, 0x3f5529b86357a1e1, 0x3bf9bb69c03848cf),
    (0x4004ce44ea1f5f32, 0x3f2ed092abd89f39, 0x3baa8d1c7b7b008b),
    (0x401399508533a95e, 0x3d929b5b0e59e22b, 0xba3fbfb12b72bd1b),
    (0x40100db4ec55041e, 0x3e4da5a3ddffe9e4, 0x3acfe056be9eeda6),
    (0x40080888bebe41fa, 0x3ef69009a4a2a7c1, 0xbb489b213dc8308f),
    (0x3ff65da40ef582e3, 0x3fa89ac4c0294785, 0xbc3b0662916500b1),
    (0x4009fb4d8a6188fc, 0x3ed2547e9d5b462c, 0x3b6e52f01808e34a),
    (0x400af4a2398615a2, 0x3ebfaa4c5491df1b, 0x3b58219a9a363175),
    (0x40064c9ba1e39094, 0x3f152f02405a85bd, 0x3bbecf4c3b101429),
    (0x4007a18e0fa9733e, 0x3efeea49dab4e7d7, 0xbb82ef504f795de7),
    (0x400d8095c0eab734, 0x3e88a123767391cb, 0xbb04e4ffd92f3b9e),
    (0x4012e97c9ab86758, 0x3db925218715debb, 0x39c0f5aac1a70278),
    (0x4005a4deee6d7cbc, 0x3f210f6073ab71bf, 0xbbc8c321a4fa8095),
    (0x3fe23e2449dc16e5, 0x3fdae3274cf05532, 0x3c70839037e7ad61),
    (0x3fe409e414339ef4, 0x3fd80dbab5fa6fc6, 0x3c0ad4073f1b4b4f),
    (0x4000f85916c89737, 0x3f661f5f650d568c, 0xbbed8306cf068afb),
    (0x3fe260b949a179b0, 0x3fdaaae69af23f6f, 0x3c769cf8bad1f8c7),
    (0x3ffb68cd57f398e3, 0x3f8f8da8964d40fb, 0xbc23257d8cd14efe),
    (0x40123394e7206430, 0x3de0f2a63b80e25a, 0xba6dee5ae9384bb1),
    (0x3ff715f0ad6f6097, 0x3fa525407ad5f519, 0x3c477771a506a703),
    (0x400b00e9be8e1ea0, 0x3ebe5c31901f30db, 0x3b5a135d91411caa),
    (0x3ff1cb8787b8f7d1, 0x3fbda1a98a0c05c7, 0x3c2ef329ed5fa1e5),
    (0x400bb4dd4d9ef2a8, 0x3eb0412451f5cb0b, 0xbb5fffb649c621e7),
    (0x3ff92074a600b7c8, 0x3f9afd437bbaf66f, 0x3c36be508df4573d),
    (0x4004648b658c616b, 0x3f3475cc08b7f908, 0xbbd71295e6b57c59),
    (0x3fd9a5d1b81f6636, 0x3fe244bca820ead4, 0x3c88c26e83fd792d),
    (0x400dda571aeb4241, 0x3e8198b62ec4ac6d, 0x3b256bbc8f401b3c),
    (0x4006a55136f6f579, 0x3f1060f0b9fbab95, 0x3bba85ea9e1ff01d),
    (0x400d4a6bb6ba8340, 0x3e8e1da664a5666f, 0xbb296fbbc77f2668),
    (0x3fe8bd56d1edf64a, 0x3fd18d2d20a223a0, 0x3c5b3e0653022319),
    (0x3ffcb11703ec0e8b, 0x3f86f6625ca7d8e4, 0xbc0ee708dd381a1e),
    (0x400abd59f02e7ff0, 0x3ec31e8bc0c289da, 0x3b6f843680a01507),
    (0x3ff9c195469d44d0, 0x3f975c5fda8bc1e3, 0x3c0b6b7cfa30be03),
    (0x4011bed4680b445b, 0x3df83038a798d8d4, 0xba7466fe037b72a7),
    (0x3f7517098a1d0780, 0x3fefd067c90b883d, 0x3c8105fd6f167d32),
    (0x3ff715123d440b0a, 0x3fa5292a5d670171, 0x3c48e9fc341b8914),
    (0x3fd29af3a8467708, 0x3fe5ca9c4278550f, 0x3c6d85214b1b094c),
    (0x40131a14a3cfbe38, 0x3dafbaa4c14b4745, 0x3a2ba071c3080178),
    (0x4005709aee8509ac, 0x3f23bd7e7565c376, 0x3bc82566570ed167),
    (0x400273e3418593bd, 0x3f52200ac328166f, 0x3bb147ddf45f34d5),
    (0x3ffe014d0fab6ec9, 0x3f8061eefdfd16f2, 0xbbd5313c55af40ea),
    (0x3ff3e94bbdca5608, 0x3fb41355a21fac7f, 0xbc577992cb937746),
    (0x400701d84a461756, 0x3f08f3f7da63cb6d, 0xbba61d5e63b471be),
    (0x401396ec5ed80f3e, 0x3d930bd2b09b4366, 0xba2c4ba778b58fd9),
    (0x400e5a969b85c344, 0x3e75a07e2aeb145f, 0x3b1c03ee6e288e5e),
    (0x400bf7c34db889fc, 0x3ea9ac5c4065dc4e, 0xbb489e10958c41c7),
    (0x40091e95b8a2651c, 0x3ee2d1967bc7b1fc, 0xbb8928f5c77b2b3d),
    (0x4007dba9d59afa50, 0x3ef9e59b78026e4a, 0xbb7351aa9c0b154b),
    (0x40124689068ce084, 0x3ddc852af69bf2e0, 0x3a7e1defc6f952df),
    (0x401133dc90190c1a, 0x3e14631ccaceae1e, 0xbaab21416bfc03a4),
    (0x40015294d314b530, 0x3f61ff58068b8a47, 0xbc004d37e1a200e2),
    (0x3fdd50e01e77ad8c, 0x3fe08c382ad4afe5, 0x3c6c67401b0202d9),
    (0x3ff03adbca45448c, 0x3fc361b6ee795042, 0x3c6f9067bb4678e5),
    (0x3fb313a32e2abc80, 0x3fed5072d4249164, 0x3c5907af24d37b0d),
    (0x4002613311ce5878, 0x3f52f7c6e98b3393, 0xbbfee0da8b5eca54),
    (0x3fec54ed978e56f2, 0x3fcaf2c79720b4d2, 0x3c6ff54267ac0145),
    (0x3ff46632a24000f4, 0x3fb245f53d895435, 0x3c22101345d612ed),
    (0x40021b935146ddab, 0x3f566fb80b461dbe, 0x3be433817a3aaf03),
    (0x400ae5224d229478, 0x3ec0b1e12f989952, 0x3b6f7382b483f766),
    (0x400023d2ec356d11, 0x3f71bae25cd5101c, 0xbbe7225cf0b6090a),
    (0x3fa3c0f873f5c5e0, 0x3fee9b8a65c8112e, 0x3c75a3c9679c9af8),
    (0x40117f3cf54ffcea, 0x3e052f8643b9cad5, 0xbaa38d334b245a4c),
    (0x4003aa6122e9f478, 0x3f40a6e0faeebf98, 0x3bd8143bc3f158c0),
    (0x40105518788df8c6, 0x3e4096388fa7e205, 0xbacf75efe04101a6),
    (0x3ff474a6e5067722, 0x3fb212d4331ed77d, 0x3c57c2f7a868ab3e),
    (0x4003ddcae2a4a36d, 0x3f3d281741538209, 0x3bdb0c7111ad0147),
    (0x40103696bf4c418b, 0x3e45489f752a204c, 0x3ae551d8278b3aa3),
    (0x400df0d071852018, 0x3e802a901ef57719, 0xbb21d5b3fb4b0dda),
    (0x400aa957b6d4b0d0, 0x3ec4768ecf26da50, 0xbb245ba46b031945),
    (0x4010f41db4605786, 0x3e21967261ab5bf0, 0x3a95e079ef708cc7),
    (0x400327acdf68197f, 0x3f4739cdd2e69be0, 0x3bd9634794468490),
    (0x3fe87dec0f55d9c1, 0x3fd1dc5f610270c9, 0x3c68d921ec927787),
    (0x3fed7e5017effd50, 0x3fc8a16764285c20, 0xbc617bbd955a9cdd),
    (0x4013fda20a1ef4b6, 0x3d7bb124572aa98d, 0x39fb06458e7f1f9b),
    (0x3fe45a5bfff9826e, 0x3fd793cb8b472c20, 0x3c782bf8b39fa818),
    (0x4003e0a1f9652503, 0x3f3cf13303b57140, 0xbbc64f381c7d42c8),
    (0x3ff22cbe88b26a44, 0x3fbbb18528729244, 0xbc49b24819b3a0eb),
    (0x40118860d41be92a, 0x3e038dfbf2c0d75d, 0xba9d09e3e2bfbff8),
    (0x3ff37c83aed31332, 0x3fb5c2bba342761a, 0x3c4529005d1b1822),
    (0x3fe5ef37e5437cf0, 0x3fd5455b8631c545, 0xbc6383a35865c5b7),
    (0x40103420694ffc5c, 0x3e45b6f68fa96ed9, 0x3ae1b1bac49227c1),
    (0x400ea24489f64897, 0x3e706b71cc874b4e, 0xbb196e709cf3d683),
    (0x3fef2dc2e731b135, 0x3fc588806558d060, 0xbc691964844df394),
    (0x4003cc2b552e8570, 0x3f3e8543cb9ba7fe, 0x3bd6300b1c38f7c2),
    (0x3ffbaea3dee423da, 0x3f8d84f6c83b8ff7, 0xbc1d5fed06d2f082),
    (0x40131094aebee91a, 0x3db15d8386cc1a91, 0x3a57796f828ed4bd),
    (0x3fef433f4a2e1e36, 0x3fc563113c2f5bd4, 0xbc654f66ac9d4805),
    (0x3fe9c86865fa8c7d, 0x3fd049f8e6d4d694, 0xbc7e77caa9e18e7a),
    (0x3fa9760453b1e8b0, 0x3fee34b437e13fd6, 0x3c89f4f60629f3ae),
    (0x4012696ae87d20e2, 0x3dd4b822db199804, 0xba531ab760202261),
    (0x40120f1a88195211, 0x3de7959a454ae929, 0xba8aad42a3e3fe0f),
    (0x3ff08ac38f63419c, 0x3fc264ff4850ff59, 0x3c6ba4db0266ede9),
    (0x400e778e5a51afb3, 0x3e735aa7a0babf52, 0xbaefea9e3f16f6f9),
    (0x4013781e2f18b461, 0x3d99b394685fdb96, 0x3a30f74e7f00debd),
    (0x40035dee92a366cb, 0x3f443f8fd694effa, 0x3bee18392651e9be),
    (0x3ff85d6a447d3290, 0x3fa003190bfa80be, 0x3c40b63145294496),
    (0x40082fe3b106d704, 0x3ef3fa7a4c048741, 0x3b8691f5688b3378),
    (0x3ff0dec760b586ca, 0x3fc1660d95a45554, 0xbc6f0ceb0279ad2e),
    (0x400ef301d29cee89, 0x3e6802c6dd7d5d3d, 0xbb037e283438d1f5),
];

fn splitmix(seed: u64, i: u64) -> u64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Uniform in [0, 1), from a seed and an index.
fn uniform(seed: u64, i: u64) -> f64 {
    (splitmix(seed, i) >> 11) as f64 / (1u64 << 53) as f64
}

/// The spacing of doubles at `t`, a normal double.
fn ulp(t: f64) -> f64 {
    let e = ((t.abs().to_bits() >> 52) as i32) - 1023;
    2f64.powi(e - 52)
}

/// The error of `v` against the exact value `hi + lo`, in ulps of the exact value.
fn ulps(v: f64, hi: f64, lo: f64) -> f64 {
    ((v - hi) - lo).abs() / ulp(hi)
}

#[test]
fn erfc_is_within_its_bound_of_mpmath() {
    let mut worst = (0.0f64, 0.0f64);
    for &(x, hi, lo) in REFERENCE.iter() {
        let (x, hi, lo) = (f64::from_bits(x), f64::from_bits(hi), f64::from_bits(lo));
        let e = ulps(erfc(x), hi, lo);
        assert!(
            e <= ERFC_ULP_BOUND,
            "erfc({x:e}) = {:e} is {e:.3} ulp from mpmath's {hi:e}",
            erfc(x)
        );
        if e > worst.0 {
            worst = (e, x);
        }
    }
    println!(
        "erfc over {} mpmath points: worst {:.4} ulp at x = {}",
        REFERENCE.len(),
        worst.0,
        worst.1
    );
    // The comparison can see an error of a few ulp: erfc moved four doubles up is outside the
    // bound at three points in four.
    let caught = REFERENCE
        .iter()
        .filter(|&&(x, hi, lo)| {
            let v = f64::from_bits(erfc(f64::from_bits(x)).to_bits() + 4);
            ulps(v, f64::from_bits(hi), f64::from_bits(lo)) > ERFC_ULP_BOUND
        })
        .count();
    assert!(
        caught * 4 > REFERENCE.len() * 3,
        "{caught} of {}",
        REFERENCE.len()
    );
}

#[test]
fn erf_and_erfc_keep_their_identities_and_special_values() {
    assert!(erfc(f64::NAN).is_nan() && erf(f64::NAN).is_nan());
    assert_eq!(erfc(f64::INFINITY), 0.0);
    assert_eq!(erfc(f64::NEG_INFINITY), 2.0);
    assert_eq!(erfc(0.0), 1.0);
    assert_eq!(erf(0.0), 0.0);
    assert_eq!(erfc(27.4), 0.0);
    assert!(erfc(26.5) > 0.0, "subnormal, not zero");
    for i in 0..2000 {
        let x = -6.0 + 12.0 * uniform(0xE1, i);
        // erf is odd; erfc(−x) = 2 − erfc(x); erf + erfc = 1. Each to the rounding of the values
        // involved, at most 4 ulp of 2 plus the bound.
        assert!((erf(-x) + erf(x)).abs() <= 4.0 * EPS, "erf odd at {x}");
        let s = erfc(x) + erfc(-x);
        assert!(
            (s - 2.0).abs() <= 8.0 * EPS,
            "erfc(x) + erfc(−x) = {s} at {x}"
        );
        assert!(
            (erf(x) + erfc(x) - 1.0).abs() <= 4.0 * EPS,
            "erf + erfc at {x}"
        );
    }
    // Across every seam of the method, on consecutive doubles, it never rises by more than the
    // two sides' bounds together: two approximations meet there, each within its bound, and
    // neither is correctly rounded, so a step of an ulp or two upwards is allowed, and counted.
    let mut rises = 0;
    for seam in [0.4375, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0] {
        let mut x = f64::from_bits(f64::to_bits(seam) - 64);
        let mut last = erfc(x);
        for _ in 0..128 {
            x = f64::from_bits(x.to_bits() + 1);
            let v = erfc(x);
            if v > last {
                rises += 1;
                assert!(
                    v - last <= 2.0 * ERFC_ULP_BOUND * ulp(v),
                    "erfc rises at {x:e} by the seam {seam}, by {:.1} ulp",
                    (v - last) / ulp(v)
                );
            }
            last = v;
        }
    }
    println!("erfc rose at {rises} of 1280 steps between consecutive doubles about the seams");
}

/// Writes `x_bits value_bits` for a million points, for `generate.py --measure`.
#[test]
#[ignore = "a dump for tools/erfc-reference/generate.py --measure, not a check"]
fn erfc_dump_for_mpmath() {
    for i in 0..1_000_000u64 {
        let x = if i % 2 == 0 {
            -3.0 + 30.0 * uniform(0xD0, i)
        } else {
            6.0 * uniform(0xD1, i)
        };
        println!("{:016x} {:016x}", x.to_bits(), erfc(x).to_bits());
    }
}

/// Charges and positions (metres) of the crystal `basis` (fractions of the cubic cell `a`, with
/// their charges) repeated `m[k]` times along axis `k`, and the box.
fn crystal(
    basis: &[([f64; 3], f64)],
    a: f64,
    m: [usize; 3],
) -> (Vec<f64>, Vec<[f64; 3]>, PeriodicBox) {
    let mut q = Vec::new();
    let mut at = Vec::new();
    for i in 0..m[0] {
        for j in 0..m[1] {
            for k in 0..m[2] {
                for &(f, c) in basis {
                    q.push(c);
                    at.push([
                        (i as f64 + f[0]) * a,
                        (j as f64 + f[1]) * a,
                        (k as f64 + f[2]) * a,
                    ]);
                }
            }
        }
    }
    (q, at, PeriodicBox::new(m.map(|n| n as f64 * a)))
}

/// Everything the Ewald sum with `p` leaves out, by brute force, joules: the real-space terms
/// `q_i q_j erfc(αr)/r` of every image pair at `r_c ≤ r < r_c + 7/α` (half of the double sum), and
/// the reciprocal terms from `k_c` to `3 k_c`, each axis with its own length.
fn omitted(q: &[f64], at: &[[f64; 3]], cell: &PeriodicBox, p: &EwaldParameters) -> f64 {
    let l = cell.lengths();
    let reach = p.cutoff + 7.0 / p.alpha;
    let images = l.map(|x| (reach / x).ceil() as i32 + 1);
    let mut real = 0.0;
    for i in 0..q.len() {
        for j in 0..q.len() {
            for nx in -images[0]..=images[0] {
                for ny in -images[1]..=images[1] {
                    for nz in -images[2]..=images[2] {
                        let n = [nx, ny, nz];
                        let d = [0, 1, 2].map(|a| at[i][a] - at[j][a] + f64::from(n[a]) * l[a]);
                        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        if r >= p.cutoff && r < reach {
                            real += 0.5 * q[i] * q[j] * erfc(p.alpha * r) / r;
                        }
                    }
                }
            }
        }
    }
    let tau = l.map(|x| 2.0 * std::f64::consts::PI / x);
    let m = tau.map(|t| (3.0 * p.k_cutoff / t).ceil() as i32);
    let a2 = 4.0 * p.alpha * p.alpha;
    let mut recip = 0.0;
    for h in -m[0]..=m[0] {
        for k in -m[1]..=m[1] {
            for n in -m[2]..=m[2] {
                let kv = [
                    tau[0] * f64::from(h),
                    tau[1] * f64::from(k),
                    tau[2] * f64::from(n),
                ];
                let k2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                if k2 <= p.k_cutoff * p.k_cutoff || k2 > 9.0 * p.k_cutoff * p.k_cutoff {
                    continue;
                }
                let (mut sr, mut si) = (0.0, 0.0);
                for (qi, r) in q.iter().zip(at) {
                    let phase = kv[0] * r[0] + kv[1] * r[1] + kv[2] * r[2];
                    sr += qi * phase.cos();
                    si += qi * phase.sin();
                }
                recip += (-k2 / a2).exp() / k2 * (sr * sr + si * si);
            }
        }
    }
    COULOMB * (real + 2.0 * std::f64::consts::PI / cell.volume() * recip)
}

fn rock_salt() -> Vec<([f64; 3], f64)> {
    let fcc = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut b: Vec<([f64; 3], f64)> = fcc.iter().map(|&f| (f, 1.0)).collect();
    b.extend(fcc.iter().map(|&f| ([f[0] + 0.5, f[1], f[2]], -1.0)));
    b
}

/// NaCl's Madelung constant on the nearest-neighbour distance, `1.74756459463318…`, and CsCl's,
/// `1.76267477307099…`, also on the nearest-neighbour distance (`a√3/2`): the energy per ion pair
/// is `−M k_e / r₀`. Each from a periodic cell at five accuracies — and rock salt again from the
/// non-cubic supercells 2×2×3 and 3×2×2, whose constant is the same, so that every axis's own
/// length is in the sum: a box that is not a cube is the only one that can tell `L_x` from `L_z`.
///
/// **Kolafa and Perram's estimate does not hold for a crystal**, and is not asserted: it assumes
/// uncorrelated charges, and an ordered lattice puts its `|S(k)|²` in Bragg peaks and its pairs in
/// shells. Measured, the error is 1.3 to 19 times the estimate on rock salt and 0.2 to 17 on
/// caesium chloride. What is asserted is exact instead: every term the two cutoffs leave out —
/// the real-space pairs beyond `r_c` over every image out to `r_c + 7/α` (`erfc 7 = 4e-23`), and
/// the wave vectors from `k_c` to `3 k_c` — summed here by brute force and added back, which
/// gives the closed form to `1e-13` at every accuracy; and the error itself falls with `δ`, to
/// under `2e-9` at the tightest (1.2e-9 for the 2×2×3 supercell, whose `k_c` is set from `V^⅓`).
#[test]
fn the_madelung_constants_of_rock_salt_and_caesium_chloride() {
    let r0 = 2.8 * ANGSTROM;
    const NACL: f64 = 1.747_564_594_633_182;
    let cs_cl = vec![([0.0, 0.0, 0.0], 1.0), ([0.5, 0.5, 0.5], -1.0)];
    let cases = [
        ("NaCl", rock_salt(), 2.0 * r0, [2usize; 3], NACL, 4.0),
        (
            "CsCl",
            cs_cl,
            2.0 * r0 / 3f64.sqrt(),
            [3; 3],
            1.762_674_773_070_988,
            1.0,
        ),
        ("NaCl 2×2×3", rock_salt(), 2.0 * r0, [2, 2, 3], NACL, 4.0),
        ("NaCl 3×2×2", rock_salt(), 2.0 * r0, [3, 2, 2], NACL, 4.0),
    ];
    for (name, basis, a, m, exact, pairs_per_cell) in cases {
        let (q, at, cell) = crystal(&basis, a, m);
        let pairs = pairs_per_cell * (m[0] * m[1] * m[2]) as f64;
        let sum_q2: f64 = q.iter().map(|x| x * x).sum();
        let cutoff = cell.half_shortest_edge();
        let mut last = f64::INFINITY;
        for accuracy in [1e-3, 1e-5, 1e-7, 1e-9, 1e-11] {
            let p = EwaldParameters::for_accuracy(&cell, cutoff, accuracy);
            let ewald = Ewald::new(cell, p);
            let e = ewald.evaluate(&q, &at, &[]).energy;
            let madelung = -e.total * r0 / (pairs * COULOMB);
            let error = (madelung - exact).abs();
            let stated = p.energy_error(&cell, sum_q2) * r0 / (pairs * COULOMB);
            let completed = -(e.total + omitted(&q, &at, &cell, &p)) * r0 / (pairs * COULOMB);
            println!(
                "{name}: δ = {accuracy:e}, α = {:.4} /Å, {} wave vectors: M = {madelung:.15}, \
                 error {error:.2e} ({:.1} times the estimate); with the omitted terms back \
                 {:.1e}",
                p.alpha * ANGSTROM,
                ewald.wave_vectors(),
                error / stated,
                completed - exact
            );
            assert!(
                (completed - exact).abs() < 1e-13,
                "{name} at δ = {accuracy:e}: completed, M = {completed} against {exact}"
            );
            assert!(
                error < last,
                "{name}: δ = {accuracy:e} is no better than the last"
            );
            last = error;
        }
        assert!(last < 2e-9, "{name}: the tightest accuracy leaves {last:e}");
    }
}

/// One charge in a cube of side `L`, neutralised by the background, is the simple cubic Wigner
/// lattice: `E = k_e q² ξ / (2L)`, `ξ = −2.837297479480620`. The background term is what makes it
/// so, at every α. A lone charge's `|S(k)|²` is `q²` at every `k`, so what the reciprocal cut
/// leaves out is exactly `(2π/V) Σ_{|k| > k_c} exp(−k²/4α²)/k²`, summed here over the wave
/// vectors out to three times `k_c`: with it added back the energy is the lattice's to rounding,
/// and the truncation is what [`EwaldParameters::reciprocal_bias`] estimates.
#[test]
fn a_lone_charge_is_the_wigner_lattice_at_every_alpha() {
    const XI: f64 = -2.837_297_479_480_62;
    let side = 15.0 * ANGSTROM;
    let cell = PeriodicBox::cubic(side);
    let exact = COULOMB * XI / (2.0 * side);
    let tau = 2.0 * std::f64::consts::PI / side;
    for cutoff_angstrom in [7.5, 6.0, 4.5] {
        for accuracy in [1e-6, 1e-9] {
            let p = EwaldParameters::for_accuracy(&cell, cutoff_angstrom * ANGSTROM, accuracy);
            let ev = Ewald::new(cell, p).evaluate(&[1.0], &[[0.3e-10, 1.1e-10, 7.0e-10]], &[]);
            let e = ev.energy;
            let m = (3.0 * p.k_cutoff / tau).ceil() as i32;
            let a2 = 4.0 * p.alpha * p.alpha;
            let mut omitted = 0.0;
            for h in -m..=m {
                for k in -m..=m {
                    for l in -m..=m {
                        let k2 = tau * tau * f64::from(h * h + k * k + l * l);
                        if k2 > p.k_cutoff * p.k_cutoff {
                            omitted += (-k2 / a2).exp() / k2;
                        }
                    }
                }
            }
            omitted *= COULOMB * 2.0 * std::f64::consts::PI / cell.volume();
            let left = e.total + omitted - exact;
            let ratio = p.reciprocal_bias(1.0) / omitted;
            println!(
                "lone charge, r_c = {cutoff_angstrom} Å, δ = {accuracy:e}: E/exact − 1 = {:.3e}; \
                 the omitted shell {:.3e}, the bias estimate {ratio:.3} of it; with the shell \
                 back, {:.1e}; background {:.4} of the total",
                e.total / exact - 1.0,
                omitted / exact.abs(),
                left / exact.abs(),
                e.background / e.total
            );
            // The nearest image is L ≥ 2 r_c away, where erfc(α L) is below 1e-30: nothing but
            // rounding is left, of the parts, whose largest is the self term.
            assert!(
                left.abs() <= 1e-13 * e.self_energy.abs(),
                "r_c {cutoff_angstrom}, δ {accuracy:e}: {left:e} left of {exact:e}"
            );
            // And the estimate is the truncation to within a factor of two either way.
            assert!(
                (0.5..2.0).contains(&ratio),
                "the bias estimate is {ratio} of the shell"
            );
            // The background is in the agreement: without it the energy would be off by a
            // million times what is left.
            assert!(e.background < 0.0 && e.background.abs() > 1e6 * left.abs());
            // No force on a lone charge.
            assert!(ev.forces[0]
                .iter()
                .all(|f| f.abs() < 1e-6 * exact.abs() / side));
        }
    }
}

/// A disordered neutral box of `n` unit-magnitude charges, at least 1.2 Å apart.
fn disordered(n: usize, side: f64, seed: u64) -> (Vec<f64>, Vec<[f64; 3]>, PeriodicBox) {
    disordered_in(n, [side; 3], seed)
}

/// The same in a box of edges `lengths`.
fn disordered_in(n: usize, lengths: [f64; 3], seed: u64) -> (Vec<f64>, Vec<[f64; 3]>, PeriodicBox) {
    let cell = PeriodicBox::new(lengths);
    let mut at: Vec<[f64; 3]> = Vec::new();
    let mut i = 0u64;
    while at.len() < n {
        let p = [0, 1, 2].map(|a| lengths[a as usize] * uniform(seed, 3 * i + a));
        i += 1;
        let close = at.iter().any(|o| {
            let d = cell.minimum_image([p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() < 1.2 * ANGSTROM
        });
        if !close {
            at.push(p);
        }
    }
    let q = (0..n)
        .map(|k| if k % 2 == 0 { 0.8 } else { -0.8 })
        .collect();
    (q, at, cell)
}

/// The total does not depend on α. At a fixed cutoff, over the accuracies, and at a fixed
/// accuracy over the cutoff — which moves α from 0.23 to 0.99 Å⁻¹ — each total is the tightest's
/// less the reciprocal bias, the signed claim: `|(E_ref − E) − (b − b_ref)| ≤ 4σ`, with `b` each
/// sum's [`EwaldParameters::reciprocal_bias`] and `σ` its [`EwaldParameters::energy_error`].
#[test]
fn the_energy_does_not_depend_on_alpha() {
    let (q, at, cell) = disordered(96, 18.0 * ANGSTROM, 0xA1FA);
    let sum_q2: f64 = q.iter().map(|x| x * x).sum();
    let (reference, reference_bias) = {
        let p = EwaldParameters::for_accuracy(&cell, 9.0 * ANGSTROM, 1e-13);
        (
            Ewald::new(cell, p).evaluate(&q, &at, &[]).energy.total,
            p.reciprocal_bias(sum_q2),
        )
    };
    let mut alphas = Vec::new();
    for (cutoff, accuracy) in [
        (9.0, 1e-3),
        (9.0, 1e-4),
        (9.0, 1e-6),
        (9.0, 1e-8),
        (9.0, 1e-10),
        (4.5, 1e-8),
        (6.0, 1e-8),
        (7.5, 1e-8),
        (4.5, 1e-11),
    ] {
        let p = EwaldParameters::for_accuracy(&cell, cutoff * ANGSTROM, accuracy);
        let e = Ewald::new(cell, p).evaluate(&q, &at, &[]).energy;
        let sigma = p.energy_error(&cell, sum_q2);
        let bias = p.reciprocal_bias(sum_q2) - reference_bias;
        let low = reference - e.total;
        println!(
            "r_c = {cutoff} Å, δ = {accuracy:e}, α = {:.4} /Å: E_ref − E = {low:.3e} J, bias \
             {bias:.3e}, left {:.2} σ (σ = {sigma:.2e})",
            p.alpha * ANGSTROM,
            (low - bias) / sigma
        );
        assert!(
            (low - bias).abs() <= 4.0 * sigma,
            "δ = {accuracy:e}: low by {low:e}, against the bias {bias:e} ± 4 × {sigma:e}"
        );
        alphas.push(p.alpha);
    }
    let (lo, hi) = alphas
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, h), &a| (l.min(a), h.max(a)));
    assert!(hi / lo > 2.0, "α spans {lo:e} to {hi:e} only");
}

/// A cluster of point charges, neutral, with dipole `μ`, in cubic boxes of growing side: the
/// periodic energy is the cluster's own Coulomb energy less `2π k_e μ²/(3V)`, the tinfoil
/// boundary's, and what is left falls as `L⁻⁵`.
#[test]
fn a_neutral_cluster_in_a_large_box_is_its_coulomb_energy_at_the_stated_rate() {
    // Six charges within 2 Å of the origin, summing to zero, with a dipole and a quadrupole.
    let q = [0.7, -0.4, -0.5, 0.3, 0.45, -0.55];
    let base = [
        [0.0, 0.0, 0.0],
        [1.2, 0.3, -0.2],
        [-0.4, 1.1, 0.5],
        [0.6, -0.9, 1.0],
        [-1.0, -0.5, -0.7],
        [0.2, 0.8, -1.3],
    ]
    .map(|p| p.map(|x| x * ANGSTROM));
    let mut plain = 0.0;
    for i in 0..6 {
        for j in i + 1..6 {
            let d = [0, 1, 2].map(|a| base[i][a] - base[j][a]);
            plain += COULOMB * q[i] * q[j] / (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        }
    }
    let mu = [0, 1, 2].map(|a| (0..6).map(|i| q[i] * base[i][a]).sum::<f64>());
    let mu2 = mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2];
    let mut rows = Vec::new();
    for side_angstrom in [12.0, 16.0, 24.0, 32.0] {
        let side = side_angstrom * ANGSTROM;
        let cell = PeriodicBox::cubic(side);
        let at: Vec<[f64; 3]> = base.iter().map(|p| p.map(|x| x + 0.37 * side)).collect();
        let p = EwaldParameters::for_accuracy(&cell, 0.5 * side, 1e-12);
        let e = Ewald::new(cell, p).evaluate(&q, &at, &[]).energy.total;
        let dipole = 2.0 * std::f64::consts::PI * COULOMB * mu2 / (3.0 * cell.volume());
        let left = e + dipole - plain;
        println!(
            "L = {side_angstrom} Å: E − E_plain = {:.4e} J, −2πμ²/3V = {:.4e}, left {left:.3e}, \
             left·L⁵ = {:.4e}",
            e - plain,
            -dipole,
            left * side_angstrom.powi(5)
        );
        rows.push((side_angstrom, e - plain, dipole, left));
    }
    for &(_, diff, dipole, left) in &rows {
        // The dipole term is most of the difference: what is left is under a tenth of it.
        assert!(left.abs() < 0.1 * dipole, "{diff:e} against {dipole:e}");
    }
    // Left × L⁵ is a constant, to the next order, which is relatively (a/L)²: measured within 2.8%
    // from 16 Å to 32 Å (5.7% at 12 Å, not compared), and held to 5%.
    let scaled: Vec<f64> = rows.iter().map(|r| r.3 * r.0.powi(5)).collect();
    for s in &scaled[1..] {
        assert!((s / scaled[3] - 1.0).abs() < 0.05, "{scaled:?}");
    }
    // And not L⁻³ or L⁻⁴: the ratio between 16 and 32 Å is 32, not 8 or 16.
    let ratio = rows[1].3 / rows[3].3;
    assert!(ratio > 24.0 && ratio < 40.0, "ratio {ratio}");
}

/// Excluded pairs are left out exactly, in any box and at any cutoff: the energy with the
/// exclusions less the energy without them is `−k_e Σ q_i q_j / r_ij` over the excluded pairs —
/// their whole Coulomb pair, which the correction's `erf(αr)/r` and the real-space `erfc(αr)/r`
/// it leaves out make together — except that a pair at or beyond `r_c` was never in the real
/// space, so for it the difference is `erf(αr)/r` alone. Held to rounding in a 40 Å box, where
/// every pair is inside `r_c`, and in an 8 Å box with `r_c` = 1.5 Å, where two are not. And in the
/// large box, beside it, the periodic energy plus the tinfoil term is the non-periodic energy
/// without the excluded pairs, to the `L⁻⁵` remainder, printed.
#[test]
fn excluded_pairs_are_left_out_exactly() {
    let q = [0.5, -0.3, -0.4, 0.2, 0.35, -0.35];
    let base = [
        [0.0, 0.0, 0.0],
        [1.0, 0.1, 0.0],
        [1.5, 1.0, 0.2],
        [-0.8, 0.6, 0.4],
        [0.3, -1.1, -0.6],
        [2.4, 1.3, 0.9],
    ]
    .map(|p| p.map(|x| x * ANGSTROM));
    let excluded = [[0, 1], [1, 2], [0, 2], [0, 3], [2, 5], [1, 5]];
    let separation = |i: usize, j: usize| {
        let d = [0, 1, 2].map(|a| base[i][a] - base[j][a]);
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    for (side, cutoff) in [(40.0, 20.0), (8.0, 1.5)] {
        let cell = PeriodicBox::cubic(side * ANGSTROM);
        let at: Vec<[f64; 3]> = base
            .iter()
            .map(|p| p.map(|x| x + 0.61 * side * ANGSTROM))
            .collect();
        let p = EwaldParameters::for_accuracy(&cell, cutoff * ANGSTROM, 1e-8);
        let ewald = Ewald::new(cell, p);
        let with = ewald.evaluate(&q, &at, &excluded).energy;
        let without = ewald.evaluate(&q, &at, &[]).energy;
        let mut expected = 0.0;
        let mut beyond = 0;
        for &[i, j] in &excluded {
            let r = separation(i, j);
            let qq = COULOMB * q[i] * q[j];
            if r < p.cutoff {
                expected -= qq / r;
            } else {
                expected -= qq * erf(p.alpha * r) / r;
                beyond += 1;
            }
        }
        let difference = with.total - without.total;
        let scale = magnitude(&with) + magnitude(&without);
        println!(
            "{side} Å box, r_c {cutoff} Å ({beyond} excluded pairs beyond it): E − E_unexcluded = \
             {difference:.12e} J against {expected:.12e}, off {:.1e} of the parts",
            (difference - expected).abs() / scale
        );
        assert!((difference - expected).abs() <= 16.0 * EPS * scale);
        assert_eq!(beyond, if cutoff < 2.0 { 2 } else { 0 });
        // The correction alone is every excluded pair's erf(αr)/r.
        let by_hand: f64 = excluded
            .iter()
            .map(|&[i, j]| {
                let r = separation(i, j);
                -COULOMB * q[i] * q[j] * erf(p.alpha * r) / r
            })
            .sum();
        assert!((with.excluded - by_hand).abs() <= 1e-13 * by_hand.abs());
        if side > 10.0 {
            let mut plain = 0.0;
            for i in 0..6 {
                for j in i + 1..6 {
                    if !excluded.contains(&[i, j]) {
                        plain += COULOMB * q[i] * q[j] / separation(i, j);
                    }
                }
            }
            let mu = [0, 1, 2].map(|a| (0..6).map(|i| q[i] * base[i][a]).sum::<f64>());
            let mu2 = mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2];
            let dipole = 2.0 * std::f64::consts::PI * COULOMB * mu2 / (3.0 * cell.volume());
            println!(
                "  and E + 2πμ²/3V − E_plain = {:.3e} J against the plain {plain:.4e} J",
                with.total + dipole - plain
            );
        }
    }
}

/// Central differences of the energy, with the tolerance the differences earn: `D(h)` and `D(2h)`
/// differ by three times the truncation term `h² E‴/6`, and rounding adds `ε_E / h` with `ε_E` the
/// energy's rounding, bounded by `64 ε` times the sum of the parts' magnitudes.
fn finite_difference(energy: impl Fn(f64) -> (f64, f64), h: f64) -> (f64, f64) {
    let d = |s: f64| {
        let (plus, scale) = energy(s);
        let (minus, _) = energy(-s);
        ((plus - minus) / (2.0 * s), scale)
    };
    let (d1, scale) = d(h);
    let (d2, _) = d(2.0 * h);
    let tolerance = 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h;
    (d1, tolerance)
}

fn magnitude(e: &pantometry_forcefield::EwaldEnergy) -> f64 {
    e.real.abs() + e.reciprocal.abs() + e.self_energy.abs() + e.excluded.abs() + e.background.abs()
}

#[test]
fn the_forces_are_the_negative_gradient() {
    let (mut q, at, cell) = disordered(40, 14.0 * ANGSTROM, 0xF0CE);
    q[0] += 0.3; // a net charge, so the background is in
    let excluded = [[0, 1], [2, 3], [3, 4], [2, 4], [10, 30]];
    for (cutoff, accuracy) in [(7.0, 1e-7), (4.5, 1e-5)] {
        let p = EwaldParameters::for_accuracy(&cell, cutoff * ANGSTROM, accuracy);
        let ewald = Ewald::new(cell, p);
        let ev = ewald.evaluate(&q, &at, &excluded);
        let h = 1e-5 * ANGSTROM;
        let mut worst = 0.0f64;
        for i in [0, 1, 2, 7, 10, 23, 39] {
            for a in 0..3 {
                let energy = |s: f64| {
                    let mut moved = at.clone();
                    moved[i][a] += s;
                    let e = ewald.evaluate(&q, &moved, &excluded).energy;
                    (e.total, magnitude(&e))
                };
                let (d, tolerance) = finite_difference(energy, h);
                let f = ev.forces[i][a];
                assert!(
                    (f + d).abs() <= tolerance,
                    "r_c {cutoff}: atom {i} axis {a}: force {f:e}, −dE/dx {:e}, tolerance \
                     {tolerance:e}",
                    -d
                );
                worst = worst.max((f + d).abs() / tolerance);
                // The tolerance is far under the force: the check can see it.
                assert!(
                    tolerance < 1e-5 * f.abs().max(1e-12),
                    "{tolerance:e} against {f:e}"
                );
            }
        }
        println!("forces, r_c = {cutoff} Å: worst |F + dE/dx| is {worst:.3} of its tolerance");
    }
}

#[test]
fn the_virial_is_the_strain_derivative() {
    let (mut q, at, cell) = disordered(40, 14.0 * ANGSTROM, 0x5EED);
    q[3] -= 0.25;
    let excluded = [[0, 1], [5, 6], [6, 7]];
    let p = EwaldParameters::for_accuracy(&cell, 6.5 * ANGSTROM, 1e-7);
    let ewald = Ewald::new(cell, p);
    let w = ewald.evaluate(&q, &at, &excluded).virial;
    let h = 1e-6;
    for a in 0..3 {
        let energy = |s: f64| {
            let mut f = [1.0; 3];
            f[a] += s;
            let scaled: Vec<[f64; 3]> = at
                .iter()
                .map(|r| [r[0] * f[0], r[1] * f[1], r[2] * f[2]])
                .collect();
            let e = ewald
                .with_cell(cell.scaled(f))
                .evaluate(&q, &scaled, &excluded)
                .energy;
            (e.total, magnitude(&e))
        };
        let (d, tolerance) = finite_difference(energy, h);
        println!(
            "virial W_{a}{a} = {:.6e} J, −dE/dε = {:.6e}, tolerance {tolerance:.2e}",
            w[a][a], -d
        );
        assert!((w[a][a] + d).abs() <= tolerance, "axis {a}");
        assert!(tolerance < 1e-6 * w[a][a].abs(), "tolerance {tolerance:e}");
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (w[a][b] - w[b][a]).abs() <= 1e-12 * w[0][0].abs(),
                "symmetric"
            );
        }
    }
}

/// The energy does not care where the system is, or which image each charge is given in.
#[test]
fn translation_and_images_change_nothing_but_rounding() {
    let (q, at, cell) = disordered(48, 13.0 * ANGSTROM, 0x1A1A);
    let excluded = [[0, 1], [1, 2], [8, 9]];
    let p = EwaldParameters::for_accuracy(&cell, 6.5 * ANGSTROM, 1e-8);
    let ewald = Ewald::new(cell, p);
    let base = ewald.evaluate(&q, &at, &excluded);
    let l = cell.lengths();
    let shift = [3.3 * ANGSTROM, -7.9 * ANGSTROM, 12.1 * ANGSTROM];
    let moved: Vec<[f64; 3]> = at
        .iter()
        .map(|r| [0, 1, 2].map(|a| r[a] + shift[a]))
        .collect();
    let images: Vec<[f64; 3]> = at
        .iter()
        .enumerate()
        .map(|(i, r)| {
            [0, 1, 2].map(|a| {
                let n = (splitmix(0x1A, 3 * i as u64 + a as u64) % 7) as f64 - 3.0;
                r[a] + n * l[a]
            })
        })
        .collect();
    let scale = magnitude(&base.energy);
    let fmax = base
        .forces
        .iter()
        .flatten()
        .fold(0.0f64, |m, f| m.max(f.abs()));
    for (what, other) in [("translated", moved), ("in other images", images)] {
        let ev = ewald.evaluate(&q, &other, &excluded);
        let de = (ev.energy.total - base.energy.total).abs();
        // Positions change by a rounding of up to 4 box lengths, 2e-15 relative; the energy by
        // that times the largest pair gradient over the pairs, held at 1e-11 of the parts.
        assert!(
            de <= 1e-11 * scale,
            "{what}: energy moved by {de:e} of {scale:e}"
        );
        for (f0, f1) in base.forces.iter().zip(&ev.forces) {
            for a in 0..3 {
                assert!((f0[a] - f1[a]).abs() <= 1e-9 * fmax, "{what}: force");
            }
        }
        println!("{what}: energy moved by {:.2e} of the parts", de / scale);
    }
}

/// The same bits on a second evaluation, and from a second sum built the same way.
#[test]
fn an_evaluation_is_its_own_bits() {
    let (q, at, cell) = disordered(30, 12.0 * ANGSTROM, 0xB175);
    let p = EwaldParameters::for_accuracy(&cell, 6.0 * ANGSTROM, 1e-6);
    let a = Ewald::new(cell, p).evaluate(&q, &at, &[[0, 1]]);
    let b = Ewald::new(cell, p).evaluate(&q, &at, &[[1, 0], [0, 1]]);
    assert_eq!(a, b);
}

/// Kolafa and Perram's real-space estimate is the RMS error it says it is, on disordered boxes:
/// with the reciprocal sum made exact (`k_c` at ten times α, where `exp(−k²/4α²)` is `e⁻²⁵`), the
/// error left by the real-space cutoff alone, over sixteen boxes, has an RMS within a factor of two
/// of the estimate, at three accuracies. The estimate is a variance over configurations; one box is
/// one draw from it, and sixteen make its RMS good to about 18%.
#[test]
fn the_real_space_estimate_is_the_error_on_disordered_boxes() {
    for accuracy in [1e-4, 1e-6, 1e-8] {
        let mut sum_sq = 0.0;
        let mut estimate = 0.0;
        for seed in 0..16u64 {
            let (q, at, cell) = disordered(96, 18.0 * ANGSTROM, 0xC0DE + seed);
            let sum_q2: f64 = q.iter().map(|x| x * x).sum();
            let cutoff = 6.0 * ANGSTROM;
            let p = EwaldParameters::for_accuracy(&cell, cutoff, accuracy);
            let exact_k = EwaldParameters {
                k_cutoff: 10.0 * p.alpha,
                ..p
            };
            let reference = EwaldParameters {
                cutoff: 9.0 * ANGSTROM,
                ..exact_k
            };
            let e = Ewald::new(cell, exact_k)
                .evaluate(&q, &at, &[])
                .energy
                .total;
            let r = Ewald::new(cell, reference)
                .evaluate(&q, &at, &[])
                .energy
                .total;
            // The reference's own real-space error at 9 Å is exp(−(9/6)² p²) of this one's.
            sum_sq += (e - r) * (e - r);
            estimate = COULOMB * sum_q2 / cutoff * p.real_space_error(&cell);
        }
        let rms = (sum_sq / 16.0).sqrt();
        println!(
            "δ = {accuracy:e}: the real-space error's RMS over 16 boxes is {:.3} of Kolafa and \
             Perram's estimate",
            rms / estimate
        );
        assert!(
            (0.5..2.0).contains(&(rms / estimate)),
            "{rms:e} against {estimate:e}"
        );
    }
}

/// The cell list finds every pair inside the cutoff and no other: the real-space energy of 120
/// charges in a box whose edges are not whole multiples of half the cutoff, 17.1 × 15.3 × 14.0 Å
/// with `r_c` = 5 Å (six, six and five cells across), against the same sum over every pair by
/// minimum image, written here. Each box elsewhere in these tests has `2L/r_c` whole, where a
/// cell count rounded the wrong way still comes out right.
#[test]
fn the_cell_list_is_every_pair() {
    let (q, at, cell) = disordered_in(120, [17.1, 15.3, 14.0].map(|x| x * ANGSTROM), 0xCE11);
    let p = EwaldParameters::for_accuracy(&cell, 5.0 * ANGSTROM, 1e-6);
    let e = Ewald::new(cell, p).evaluate(&q, &at, &[]).energy;
    let mut brute = 0.0;
    let mut count = 0;
    for i in 0..q.len() {
        for j in i + 1..q.len() {
            let d = cell.minimum_image([0, 1, 2].map(|a| at[i][a] - at[j][a]));
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if r < p.cutoff {
                brute += COULOMB * q[i] * q[j] * erfc(p.alpha * r) / r;
                count += 1;
            }
        }
    }
    let scale: f64 = (0..q.len())
        .flat_map(|i| (i + 1..q.len()).map(move |j| (i, j)))
        .map(|(i, j)| {
            let d = cell.minimum_image([0, 1, 2].map(|a| at[i][a] - at[j][a]));
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if r < p.cutoff {
                (COULOMB * q[i] * q[j] * erfc(p.alpha * r) / r).abs()
            } else {
                0.0
            }
        })
        .sum();
    println!(
        "cell list: {count} pairs inside r_c; real space {:.12e} J against {brute:.12e}, off {:.1e} \
         of the pairs' magnitudes",
        e.real,
        (e.real - brute).abs() / scale
    );
    assert!((e.real - brute).abs() <= 16.0 * EPS * scale);
}

/// The off-diagonal virial against a sheared box. A simple shear `x → x + ε y` takes every atom and
/// the box with it, keeps the volume and every phase `k·r`, and moves each wave vector to
/// `k_y − ε k_x`, so the energy of the sheared box is written here as a sum over the same pairs
/// (inside `r_c` before the shear, none within the step of it) and the same integer wave vectors,
/// differentiated by central differences, and set against `−W_xy`, `−W_xz` and `−W_yz`: each to
/// the tolerance the differences earn, in a box that is not a cube, with exclusions and a net
/// charge.
#[test]
fn the_off_diagonal_virial_is_the_shear_derivative() {
    let (mut q, at, cell) = disordered_in(48, [13.0, 11.5, 12.2].map(|x| x * ANGSTROM), 0x5AEA);
    q[5] += 0.2;
    let excluded = [[0usize, 1usize], [2, 3], [3, 4]];
    let p = EwaldParameters::for_accuracy(&cell, 5.5 * ANGSTROM, 1e-7);
    let ev = Ewald::new(cell, p).evaluate(&q, &at, &excluded);
    let l = cell.lengths();
    let pi = std::f64::consts::PI;
    // The pairs: minimum-image separations before the shear.
    let mut pairs = Vec::new();
    for i in 0..q.len() {
        for j in i + 1..q.len() {
            let d = cell.minimum_image([0, 1, 2].map(|a| at[i][a] - at[j][a]));
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let ex = excluded.contains(&[i, j]);
            if ex || r < p.cutoff {
                assert!(
                    ex || (r - p.cutoff).abs() > 1e-5 * p.cutoff,
                    "a pair at the cutoff"
                );
                pairs.push((q[i] * q[j], d, ex));
            }
        }
    }
    // The wave vectors, as `Ewald::new` chooses them, with |S(k)|² from the positions.
    let extent = [0, 1, 2].map(|a| (p.k_cutoff * l[a] / (2.0 * pi)).floor() as i32);
    let mut waves = Vec::new();
    for h in 0..=extent[0] {
        for k in -extent[1]..=extent[1] {
            for m in -extent[2]..=extent[2] {
                if !(h > 0 || (h == 0 && (k > 0 || (k == 0 && m > 0)))) {
                    continue;
                }
                let kv = [0, 1, 2].map(|a| 2.0 * pi * f64::from([h, k, m][a]) / l[a]);
                if kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2] > p.k_cutoff * p.k_cutoff {
                    continue;
                }
                let (mut sr, mut si) = (0.0, 0.0);
                for (qi, r) in q.iter().zip(&at) {
                    let phase = kv[0] * r[0] + kv[1] * r[1] + kv[2] * r[2];
                    sr += qi * phase.cos();
                    si += qi * phase.sin();
                }
                waves.push((kv, sr * sr + si * si));
            }
        }
    }
    // The energy of the box sheared by ε, `a` moving with `b`, less the terms a shear leaves
    // alone (self and background), and its parts' magnitudes.
    let sheared = |a: usize, b: usize, eps: f64| -> (f64, f64) {
        let (mut e, mut scale) = (0.0, 0.0);
        for &(qq, d, ex) in &pairs {
            let mut d = d;
            d[a] += eps * d[b];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let t = if ex {
                -COULOMB * qq * erf(p.alpha * r) / r
            } else {
                COULOMB * qq * erfc(p.alpha * r) / r
            };
            e += t;
            scale += t.abs();
        }
        let a2 = 4.0 * p.alpha * p.alpha;
        for &(k0, s2) in &waves {
            let mut k = k0;
            k[b] -= eps * k0[a];
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            let t = COULOMB * 4.0 * pi / cell.volume() * (-k2 / a2).exp() / k2 * s2;
            e += t;
            scale += t.abs();
        }
        (e, scale)
    };
    // At no shear the written sum is the Ewald sum's, less its self and background terms.
    let (e0, scale0) = sheared(0, 1, 0.0);
    let ewald = ev.energy.real + ev.energy.reciprocal + ev.energy.excluded;
    assert!(
        (e0 - ewald).abs() <= 1e-12 * scale0,
        "{e0:e} against {ewald:e}"
    );
    let h = 1e-6;
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let (d, tolerance) = finite_difference(|s| sheared(a, b, s), h);
        // dU/dε = Σ ∂U/∂r_a r_b = −Σ F_a r_b = −W_ba, and W is symmetric.
        let w = ev.virial[b][a];
        println!(
            "W_{a}{b} = {w:.8e} J, −dU/dε = {:.8e}, tolerance {tolerance:.1e}",
            -d
        );
        assert!((w + d).abs() <= tolerance, "W_{a}{b}");
        assert!(
            tolerance < 1e-4 * w.abs(),
            "the check can see W_{a}{b}: {tolerance:e}"
        );
    }
}
