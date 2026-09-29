pub const HBFilter = struct { taps: []const f64, atten: f64 };

const HBKernel_4A_taps = [_]f64{
    6.1729335650971517e-001,
    -1.5963945620743250e-001,
    5.5073370934086312e-002,
    -1.4603578989932850e-002,
};
pub const HBKernel_4A = HBFilter{ .taps = &HBKernel_4A_taps, .atten = 54.5176 };

const HBKernel_5A_taps = [_]f64{
    6.2068807424902472e-001,
    -1.6827573634467302e-001,
    6.5263016720721170e-002,
    -2.2483331611592005e-002,
    5.2917326684281110e-003,
};
pub const HBKernel_5A = HBFilter{ .taps = &HBKernel_5A_taps, .atten = 66.3075 };

const HBKernel_6A_taps = [_]f64{
    6.2187202340480707e-001,
    -1.7132842113816371e-001,
    6.9019169178765674e-002,
    -2.5799728312695277e-002,
    7.4880112525741666e-003,
    -1.2844465869952567e-003,
};
pub const HBKernel_6A = HBFilter{ .taps = &HBKernel_6A_taps, .atten = 89.5271 };

const HBKernel_7A_taps = [_]f64{
    6.2354494135775851e-001,
    -1.7571220703702045e-001,
    7.4529843603968457e-002,
    -3.0701736822442153e-002,
    1.0716755639039573e-002,
    -2.7833422930759735e-003,
    4.1118797093875510e-004,
};
pub const HBKernel_7A = HBFilter{ .taps = &HBKernel_7A_taps, .atten = 105.2842 };

const HBKernel_8A_taps = [_]f64{
    6.2488363107953926e-001,
    -1.7924942606514119e-001,
    7.9068155655640557e-002,
    -3.4907523415495731e-002,
    1.3710256799907897e-002,
    -4.3991142586987933e-003,
    1.0259190163889602e-003,
    -1.3278941979339359e-004,
};
pub const HBKernel_8A = HBFilter{ .taps = &HBKernel_8A_taps, .atten = 121.0063 };

const HBKernel_9A_taps = [_]f64{
    6.2597763804021977e-001,
    -1.8216414325139055e-001,
    8.2879104876726728e-002,
    -3.8563442248249404e-002,
    1.6471530499739394e-002,
    -6.0489108881335227e-003,
    1.7805283804140392e-003,
    -3.7533200112729561e-004,
    4.3172840558735476e-005,
};
pub const HBKernel_9A = HBFilter{ .taps = &HBKernel_9A_taps, .atten = 136.6982 };

const HBKernel_10A_taps = [_]f64{
    6.2688767582974092e-001,
    -1.8460766807559420e-001,
    8.6128943000481864e-002,
    -4.1774474147006607e-002,
    1.9014801985747346e-002,
    -7.6870397465866507e-003,
    2.6264590175341853e-003,
    -7.1106660285478562e-004,
    1.3645852036179345e-004,
    -1.4113888783332969e-005,
};
pub const HBKernel_10A = HBFilter{ .taps = &HBKernel_10A_taps, .atten = 152.3572 };

const HBKernel_11A_taps = [_]f64{
    6.2667167706948146e-001,
    -1.8407153342635879e-001,
    8.5529995610836046e-002,
    -4.1346831462361310e-002,
    1.8844831691322637e-002,
    -7.7125170365394992e-003,
    2.7268674860562087e-003,
    -7.9745028501057233e-004,
    1.8116344606360795e-004,
    -2.8569149754241848e-005,
    2.3667022010173616e-006,
};
pub const HBKernel_11A = HBFilter{ .taps = &HBKernel_11A_taps, .atten = 183.7962 };

const HBKernel_12A_taps = [_]f64{
    6.2747849730367999e-001,
    -1.8623616784506747e-001,
    8.8409755898467945e-002,
    -4.4207468821462342e-002,
    2.1149175945115381e-002,
    -9.2551508371115209e-003,
    3.5871562170822330e-003,
    -1.1923167653750219e-003,
    3.2627812189920129e-004,
    -6.9106902511490413e-005,
    1.0122897863125124e-005,
    -7.7531878906846174e-007,
};
pub const HBKernel_12A = HBFilter{ .taps = &HBKernel_12A_taps, .atten = 199.4768 };

const HBKernel_13A_taps = [_]f64{
    6.2816416252367324e-001,
    -1.8809076955230414e-001,
    9.0918539867353029e-002,
    -4.6765502683599310e-002,
    2.3287520498995663e-002,
    -1.0760627245014184e-002,
    4.4853922948425683e-003,
    -1.6438775426910800e-003,
    5.1441312354764978e-004,
    -1.3211725685765050e-004,
    2.6191319837779187e-005,
    -3.5802430606313093e-006,
    2.5491278270628601e-007,
};
pub const HBKernel_13A = HBFilter{ .taps = &HBKernel_13A_taps, .atten = 215.1364 };

const HBKernel_14A_taps = [_]f64{
    6.2875473120929948e-001,
    -1.8969941936903847e-001,
    9.3126094480960403e-002,
    -4.9067251179869126e-002,
    2.5273008851199916e-002,
    -1.2218646153393291e-002,
    5.4048942085580280e-003,
    -2.1409919546078581e-003,
    7.4250292812927973e-004,
    -2.1924542206832172e-004,
    5.3015808983125091e-005,
    -9.8743034923598196e-006,
    1.2650391141650221e-006,
    -8.4146674637474946e-008,
};
pub const HBKernel_14A = HBFilter{ .taps = &HBKernel_14A_taps, .atten = 230.7526 };

const HBKernel_2B_taps = [_]f64{
    5.7361525854329076e-001,
    -7.5092074924827903e-002,
};
pub const HBKernel_2B = HBFilter{ .taps = &HBKernel_2B_taps, .atten = 56.6007 };

const HBKernel_3B_taps = [_]f64{
    5.9277038608066912e-001,
    -1.0851340190268854e-001,
    1.5813570475513079e-002,
};
pub const HBKernel_3B = HBFilter{ .taps = &HBKernel_3B_taps, .atten = 83.0295 };

const HBKernel_4B_taps = [_]f64{
    6.0140277542879617e-001,
    -1.2564483854574138e-001,
    2.7446500598038322e-002,
    -3.2051079559057435e-003,
};
pub const HBKernel_4B = HBFilter{ .taps = &HBKernel_4B_taps, .atten = 123.4724 };

const HBKernel_5B_taps = [_]f64{
    6.0818642429088932e-001,
    -1.3981140187175697e-001,
    3.8489164054503623e-002,
    -7.6218861797853104e-003,
    7.5772358130952392e-004,
};
pub const HBKernel_5B = HBFilter{ .taps = &HBKernel_5B_taps, .atten = 152.4411 };

const HBKernel_6B_taps = [_]f64{
    6.1278392271464355e-001,
    -1.5000053762513338e-001,
    4.7575323511364960e-002,
    -1.2320702802243476e-002,
    2.1462442592348487e-003,
    -1.8425092381892940e-004,
};
pub const HBKernel_6B = HBFilter{ .taps = &HBKernel_6B_taps, .atten = 181.2501 };

const HBKernel_7B_taps = [_]f64{
    6.1610372263478952e-001,
    -1.5767891882524138e-001,
    5.5089691170294691e-002,
    -1.6895755656366061e-002,
    3.9416643438213977e-003,
    -6.0603623791604668e-004,
    4.5632602433393365e-005,
};
pub const HBKernel_7B = HBFilter{ .taps = &HBKernel_7B_taps, .atten = 209.9472 };

const HBKernel_8B_taps = [_]f64{
    6.1861282914465976e-001,
    -1.6367179451225150e-001,
    6.1369861342939716e-002,
    -2.1184466539006987e-002,
    5.9623357510842061e-003,
    -1.2483098507454090e-003,
    1.7099297537964702e-004,
    -1.1448313239478885e-005,
};
pub const HBKernel_8B = HBFilter{ .taps = &HBKernel_8B_taps, .atten = 238.5616 };

const HBKernel_2C_taps = [_]f64{
    5.6430278013478008e-001,
    -6.4338068855763375e-002,
};
pub const HBKernel_2C = HBFilter{ .taps = &HBKernel_2C_taps, .atten = 89.0473 };

const HBKernel_3C_taps = [_]f64{
    5.8706402915551448e-001,
    -9.9362380958670449e-002,
    1.2298637065869358e-002,
};
pub const HBKernel_3C = HBFilter{ .taps = &HBKernel_3C_taps, .atten = 130.8951 };

const HBKernel_4C_taps = [_]f64{
    5.9896586134984675e-001,
    -1.2111680603434927e-001,
    2.4763118076458895e-002,
    -2.6121758132212989e-003,
};
pub const HBKernel_4C = HBFilter{ .taps = &HBKernel_4C_taps, .atten = 172.3192 };

const HBKernel_5C_taps = [_]f64{
    6.0626808285230716e-001,
    -1.3588224032740795e-001,
    3.5544305238309003e-002,
    -6.5127022377289654e-003,
    5.8255449565950768e-004,
};
pub const HBKernel_5C = HBFilter{ .taps = &HBKernel_5C_taps, .atten = 213.4984 };

const HBKernel_6C_taps = [_]f64{
    6.1120171263351242e-001,
    -1.4654486853757870e-001,
    4.4582959299131253e-002,
    -1.0840543858123995e-002,
    1.7343706485509962e-003,
    -1.3363018567985596e-004,
};
pub const HBKernel_6C = HBFilter{ .taps = &HBKernel_6C_taps, .atten = 254.5186 };

const HBKernel_1D_taps = [_]f64{
    5.0188900022775451e-001,
};
pub const HBKernel_1D = HBFilter{ .taps = &HBKernel_1D_taps, .atten = 54.4754 };

const HBKernel_2D_taps = [_]f64{
    5.6295152180538044e-001,
    -6.2953706070191726e-002,
};
pub const HBKernel_2D = HBFilter{ .taps = &HBKernel_2D_taps, .atten = 113.2139 };

const HBKernel_3D_taps = [_]f64{
    5.8621968728755036e-001,
    -9.8080551656524531e-002,
    1.1860868761997080e-002,
};
pub const HBKernel_3D = HBFilter{ .taps = &HBKernel_3D_taps, .atten = 167.1447 };

const HBKernel_4D_taps = [_]f64{
    5.9835028657163591e-001,
    -1.1999986086623511e-001,
    2.4132530854004228e-002,
    -2.4829565686819706e-003,
};
pub const HBKernel_4D = HBFilter{ .taps = &HBKernel_4D_taps, .atten = 220.6519 };

const HBKernel_1E_taps = [_]f64{
    5.0047102586416625e-001,
};
pub const HBKernel_1E = HBFilter{ .taps = &HBKernel_1E_taps, .atten = 66.5391 };

const HBKernel_2E_taps = [_]f64{
    5.6261293163933568e-001,
    -6.2613067826620017e-002,
};
pub const HBKernel_2E = HBFilter{ .taps = &HBKernel_2E_taps, .atten = 137.3173 };

const HBKernel_3E_taps = [_]f64{
    5.8600808139396787e-001,
    -9.7762185880067784e-002,
    1.1754104554493029e-002,
};
pub const HBKernel_3E = HBFilter{ .taps = &HBKernel_3E_taps, .atten = 203.2997 };

const HBKernel_4E_taps = [_]f64{
    5.9819599352772002e-001,
    -1.1972157555011861e-001,
    2.3977305567947922e-002,
    -2.4517235455853992e-003,
};
pub const HBKernel_4E = HBFilter{ .taps = &HBKernel_4E_taps, .atten = 268.8550 };

const HBKernel_1F_taps = [_]f64{
    5.0007530666642896e-001,
};
pub const HBKernel_1F = HBFilter{ .taps = &HBKernel_1F_taps, .atten = 82.4633 };

const HBKernel_2F_taps = [_]f64{
    5.6252823610146030e-001,
    -6.2528244608044792e-002,
};
pub const HBKernel_2F = HBFilter{ .taps = &HBKernel_2F_taps, .atten = 161.4049 };

const HBKernel_3F_taps = [_]f64{
    5.8595514744674237e-001,
    -9.7682725156791952e-002,
    1.1727577711117231e-002,
};
pub const HBKernel_3F = HBFilter{ .taps = &HBKernel_3F_taps, .atten = 239.4313 };

const HBKernel_1G_taps = [_]f64{
    5.0001882524896712e-001,
};
pub const HBKernel_1G = HBFilter{ .taps = &HBKernel_1G_taps, .atten = 94.5052 };

const HBKernel_2G_taps = [_]f64{
    5.6250705922479682e-001,
    -6.2507059756378394e-002,
};
pub const HBKernel_2G = HBFilter{ .taps = &HBKernel_2G_taps, .atten = 185.4886 };

const HBKernel_3G_taps = [_]f64{
    5.8594191201187384e-001,
    -9.7662868266991207e-002,
    1.1720956255134043e-002,
};
pub const HBKernel_3G = HBFilter{ .taps = &HBKernel_3G_taps, .atten = 275.5501 };

pub const FltGroupA = [_]HBFilter{
    HBKernel_4A,
    HBKernel_5A,
    HBKernel_6A,
    HBKernel_7A,
    HBKernel_8A,
    HBKernel_9A,
    HBKernel_10A,
    HBKernel_11A,
    HBKernel_12A,
    HBKernel_13A,
    HBKernel_14A,
};

pub const FltGroupB = [_]HBFilter{
    HBKernel_2B,
    HBKernel_3B,
    HBKernel_4B,
    HBKernel_5B,
    HBKernel_6B,
    HBKernel_7B,
    HBKernel_8B,
};

pub const FltGroupC = [_]HBFilter{
    HBKernel_2C,
    HBKernel_3C,
    HBKernel_4C,
    HBKernel_5C,
    HBKernel_6C,
};

pub const FltGroupD = [_]HBFilter{
    HBKernel_1D,
    HBKernel_2D,
    HBKernel_3D,
    HBKernel_4D,
};

pub const FltGroupE = [_]HBFilter{
    HBKernel_1E,
    HBKernel_2E,
    HBKernel_3E,
    HBKernel_4E,
};

pub const FltGroupF = [_]HBFilter{
    HBKernel_1F,
    HBKernel_2F,
    HBKernel_3F,
};

pub const FltGroupG = [_]HBFilter{
    HBKernel_1G,
    HBKernel_2G,
    HBKernel_3G,
};

const HBKernel_3A_Third_taps = [_]f64{
    5.9811355069551475e-001,
    -1.1793396656733847e-001,
    2.0300557211946322e-002,
};
pub const HBKernel_3A_Third = HBFilter{ .taps = &HBKernel_3A_Third_taps, .atten = 66.3726 };

const HBKernel_4A_Third_taps = [_]f64{
    6.0645499250612578e-001,
    -1.3555496505481171e-001,
    3.4022804962365975e-002,
    -4.9535418595798757e-003,
};
pub const HBKernel_4A_Third = HBFilter{ .taps = &HBKernel_4A_Third_taps, .atten = 90.2546 };

const HBKernel_5A_Third_taps = [_]f64{
    6.1014115058940210e-001,
    -1.4393081816629907e-001,
    4.1760642892852244e-002,
    -8.9692183234056175e-003,
    9.9871340618342070e-004,
};
pub const HBKernel_5A_Third = HBFilter{ .taps = &HBKernel_5A_Third_taps, .atten = 126.5507 };

const HBKernel_6A_Third_taps = [_]f64{
    6.1439563420546972e-001,
    -1.5360187826905250e-001,
    5.0840891345687034e-002,
    -1.4053648740561121e-002,
    2.6771286587305727e-003,
    -2.5815816044823123e-004,
};
pub const HBKernel_6A_Third = HBFilter{ .taps = &HBKernel_6A_Third_taps, .atten = 150.1839 };

const HBKernel_7A_Third_taps = [_]f64{
    6.1747493476329918e-001,
    -1.6087373733313212e-001,
    5.8263075641409430e-002,
    -1.8872408173431318e-002,
    4.7421376543513687e-003,
    -8.0196529612267474e-004,
    6.7964807393798996e-005,
};
pub const HBKernel_7A_Third = HBFilter{ .taps = &HBKernel_7A_Third_taps, .atten = 173.7068 };

const HBKernel_8A_Third_taps = [_]f64{
    6.1980610947775050e-001,
    -1.6654070578314714e-001,
    6.4416567441730327e-002,
    -2.3307744348719822e-002,
    6.9909157372312443e-003,
    -1.5871946293364403e-003,
    2.4017727382382763e-004,
    -1.8125308241541697e-005,
};
pub const HBKernel_8A_Third = HBFilter{ .taps = &HBKernel_8A_Third_taps, .atten = 197.1454 };

const HBKernel_9A_Third_taps = [_]f64{
    6.2163188951899306e-001,
    -1.7108115323810941e-001,
    6.9588370095600260e-002,
    -2.7339625080613838e-002,
    9.2954469183791771e-003,
    -2.5537179959555429e-003,
    5.2572290897951021e-004,
    -7.1813356135154921e-005,
    4.8802382808892154e-006,
};
pub const HBKernel_9A_Third = HBFilter{ .taps = &HBKernel_9A_Third_taps, .atten = 220.5199 };

const HBKernel_2B_Third_taps = [_]f64{
    5.6748544264806311e-001,
    -6.7764090509431732e-002,
};
pub const HBKernel_2B_Third = HBFilter{ .taps = &HBKernel_2B_Third_taps, .atten = 71.0965 };

const HBKernel_3B_Third_taps = [_]f64{
    5.8793612182667199e-001,
    -1.0070583248877293e-001,
    1.2771337947163834e-002,
};
pub const HBKernel_3B_Third = HBFilter{ .taps = &HBKernel_3B_Third_taps, .atten = 115.7707 };

const HBKernel_4B_Third_taps = [_]f64{
    5.9960155600862808e-001,
    -1.2228154335199336e-001,
    2.5433718917694709e-002,
    -2.7537562530837154e-003,
};
pub const HBKernel_4B_Third = HBFilter{ .taps = &HBKernel_4B_Third_taps, .atten = 152.1535 };

const HBKernel_5B_Third_taps = [_]f64{
    6.0676859170554343e-001,
    -1.3689667009876413e-001,
    3.6288512631926818e-002,
    -6.7838855305035351e-003,
    6.2345167677087547e-004,
};
pub const HBKernel_5B_Third = HBFilter{ .taps = &HBKernel_5B_Third_taps, .atten = 188.2914 };

const HBKernel_6B_Third_taps = [_]f64{
    6.1161456341904397e-001,
    -1.4743901958274458e-001,
    4.5344160157313275e-002,
    -1.1207371780924531e-002,
    1.8328497112594935e-003,
    -1.4518193006359589e-004,
};
pub const HBKernel_6B_Third = HBFilter{ .taps = &HBKernel_6B_Third_taps, .atten = 224.2705 };

const HBKernel_1C_Third_taps = [_]f64{
    5.0336730531430562e-001,
};
pub const HBKernel_1C_Third = HBFilter{ .taps = &HBKernel_1C_Third_taps, .atten = 49.4544 };

const HBKernel_2C_Third_taps = [_]f64{
    5.6330232648142819e-001,
    -6.3309247177420452e-002,
};
pub const HBKernel_2C_Third = HBFilter{ .taps = &HBKernel_2C_Third_taps, .atten = 103.1970 };

const HBKernel_3C_Third_taps = [_]f64{
    5.8643891113580415e-001,
    -9.8411593011583087e-002,
    1.1972706651483846e-002,
};
pub const HBKernel_3C_Third = HBFilter{ .taps = &HBKernel_3C_Third_taps, .atten = 152.1195 };

const HBKernel_4C_Third_taps = [_]f64{
    5.9851012363917222e-001,
    -1.2028885239978220e-001,
    2.4294521083140615e-002,
    -2.5157924156609776e-003,
};
pub const HBKernel_4C_Third = HBFilter{ .taps = &HBKernel_4C_Third_taps, .atten = 200.6182 };

const HBKernel_5C_Third_taps = [_]f64{
    6.0590922882030196e-001,
    -1.3515953438018685e-001,
    3.5020857107815606e-002,
    -6.3256196990467053e-003,
    5.5506815147598793e-004,
};
pub const HBKernel_5C_Third = HBFilter{ .taps = &HBKernel_5C_Third_taps, .atten = 248.8730 };

const HBKernel_1D_Third_taps = [_]f64{
    5.0083794231068057e-001,
};
pub const HBKernel_1D_Third = HBFilter{ .taps = &HBKernel_1D_Third_taps, .atten = 61.5357 };

const HBKernel_2D_Third_taps = [_]f64{
    5.6270074379958690e-001,
    -6.2701174487726344e-002,
};
pub const HBKernel_2D_Third = HBFilter{ .taps = &HBKernel_2D_Third_taps, .atten = 127.3167 };

const HBKernel_3D_Third_taps = [_]f64{
    5.8606296210323228e-001,
    -9.7844644765123029e-002,
    1.1781683046528768e-002,
};
pub const HBKernel_3D_Third = HBFilter{ .taps = &HBKernel_3D_Third_taps, .atten = 188.2990 };

const HBKernel_4D_Third_taps = [_]f64{
    5.9823601243162516e-001,
    -1.1979368994739022e-001,
    2.4017458606412575e-002,
    -2.4597810910081913e-003,
};
pub const HBKernel_4D_Third = HBFilter{ .taps = &HBKernel_4D_Third_taps, .atten = 248.8580 };

const HBKernel_1E_Third_taps = [_]f64{
    5.0013388897382527e-001,
};
pub const HBKernel_1E_Third = HBFilter{ .taps = &HBKernel_1E_Third_taps, .atten = 77.4651 };

const HBKernel_2E_Third_taps = [_]f64{
    5.6255019604317880e-001,
    -6.2550222932381064e-002,
};
pub const HBKernel_2E_Third = HBFilter{ .taps = &HBKernel_2E_Third_taps, .atten = 151.4084 };

const HBKernel_3E_Third_taps = [_]f64{
    5.8596887234201078e-001,
    -9.7703321113080305e-002,
    1.1734448777069783e-002,
};
pub const HBKernel_3E_Third = HBFilter{ .taps = &HBKernel_3E_Third_taps, .atten = 224.4365 };

const HBKernel_1F_Third_taps = [_]f64{
    5.0003346776264190e-001,
};
pub const HBKernel_1F_Third = HBFilter{ .taps = &HBKernel_1F_Third_taps, .atten = 89.5075 };

const HBKernel_2F_Third_taps = [_]f64{
    5.6251254964097952e-001,
    -6.2512551321105267e-002,
};
pub const HBKernel_2F_Third = HBFilter{ .taps = &HBKernel_2F_Third_taps, .atten = 175.4932 };

const HBKernel_3F_Third_taps = [_]f64{
    5.8594534336747051e-001,
    -9.7668015838639821e-002,
    1.1722672471262996e-002,
};
pub const HBKernel_3F_Third = HBFilter{ .taps = &HBKernel_3F_Third_taps, .atten = 260.5645 };

const HBKernel_1G_Third_taps = [_]f64{
    5.0000836666064941e-001,
};
pub const HBKernel_1G_Third = HBFilter{ .taps = &HBKernel_1G_Third_taps, .atten = 101.5490 };

const HBKernel_2G_Third_taps = [_]f64{
    5.6250313744943459e-001,
    -6.2503137554435345e-002,
};
pub const HBKernel_2G_Third = HBFilter{ .taps = &HBKernel_2G_Third_taps, .atten = 199.5761 };

const HBKernel_3G_Third_taps = [_]f64{
    5.8593945786963764e-001,
    -9.7659186853499613e-002,
    1.1719728983863425e-002,
};
pub const HBKernel_3G_Third = HBFilter{ .taps = &HBKernel_3G_Third_taps, .atten = 296.5185 };

pub const FltGroupA_Third = [_]HBFilter{
    HBKernel_3A_Third,
    HBKernel_4A_Third,
    HBKernel_5A_Third,
    HBKernel_6A_Third,
    HBKernel_7A_Third,
    HBKernel_8A_Third,
    HBKernel_9A_Third,
};

pub const FltGroupB_Third = [_]HBFilter{
    HBKernel_2B_Third,
    HBKernel_3B_Third,
    HBKernel_4B_Third,
    HBKernel_5B_Third,
    HBKernel_6B_Third,
};

pub const FltGroupC_Third = [_]HBFilter{
    HBKernel_1C_Third,
    HBKernel_2C_Third,
    HBKernel_3C_Third,
    HBKernel_4C_Third,
    HBKernel_5C_Third,
};

pub const FltGroupD_Third = [_]HBFilter{
    HBKernel_1D_Third,
    HBKernel_2D_Third,
    HBKernel_3D_Third,
    HBKernel_4D_Third,
};

pub const FltGroupE_Third = [_]HBFilter{
    HBKernel_1E_Third,
    HBKernel_2E_Third,
    HBKernel_3E_Third,
};

pub const FltGroupF_Third = [_]HBFilter{
    HBKernel_1F_Third,
    HBKernel_2F_Third,
    HBKernel_3F_Third,
};

pub const FltGroupG_Third = [_]HBFilter{
    HBKernel_1G_Third,
    HBKernel_2G_Third,
    HBKernel_3G_Third,
};

