########################################################################################################################
#                                                                                                                      #
#                                    Software Projekt: Frechet Distanz                                                 #
#                                    Teilgebiet: Ellipsen-Alg. einer Zelle                                             #
#                                    Erstellt: WS 16/17 FU Berlin                                                      #
#                                                                                                                      #
#                                    Team: Josephine Mertens, Jana Kirschner,                                          #
#                                    Alexander Korzech, Fabian Kovacs, Alexander                                       #
#                                    Timme, Kilian Kraatz & Anton Begehr                                               #
#                                                                                                                      #
########################################################################################################################

# -*- coding: utf-8 -*-

from Graphics import *


n_l = 7  # number of height lines
n_p = 50  # points per ellipsis

'''# path_p:
p_xs = [0, 1]
p_ys = [0, 1]
# path_q:
q_xs = [1, 0]
q_ys = [0, 0]'''

'''# path_p:
p_xs = [-8, 12, 12, 0, -10]
p_ys = [-2, -2, 2, 10, 5]
# path_q:
q_xs = [-10, 10, -10]
q_ys = [0, 0, 4]'''

'''# new type of events ex. 1
# path_p:
p_xs = [0, 1.5, 1, 2]
p_ys = [0, 0, 1+math.sqrt(0.5), 1+math.sqrt(0.5)]
# path_q:
q_xs = [0, 2, 1]
q_ys = [math.sqrt(0.5), math.sqrt(0.5), 0]'''

'''# new type of event ex. 2
# path_p:
p_xs = [4, 1.7, 2, 0]
p_ys = [1-math.sqrt(0.5), 1-math.sqrt(0.5), 0, 0]
# path_q:
q_xs = [4, 0]
q_ys = [1, 1]'''

'''# old type of event acute:
# path_p:
p_xs = [0, 10]
p_ys = [0, 0]
# path_q:
q_xs = [1, 4, 7, 10]
q_ys = [0, 3, -3, 0]'''

'''# old type of event not acute:
# path_p:
p_xs = [0, 10]
p_ys = [0, 0]
# path_q:
q_xs = [1, 7, 4, 10]
q_ys = [0, 3, -3, 0]'''

'''# critical events test 1
# path_p:
p_xs = [5, 10, 10, 0, 0, 10]
p_ys = [10, 10, 0, 0, 5, 5]
# path_q:
q_xs = [-10, -5, -5, -15, -15, -5]
q_ys = [10, 10, 5, 5, 0, 0]'''

'''# critical events test 2
# path_p:
p_xs = [0, 20]
p_ys = [0, 0]
# path_q:
q_xs = [10, 15, 20]
q_ys = [0, 10, 0]'''

'''# critical events test 3
# path_p:
p_xs = [0, 10, 3, 0, 10]
p_ys = [0, 0, 8, 10, 10]
# path_q:
q_xs = [0, 10]
q_ys = [3, 3]'''

'''# random test 1
# path_p:
p_xs = np.random.random_sample(9)
p_ys = np.random.random_sample(9)
# path_q:
q_xs = np.random.random_sample(9)
q_ys = np.random.random_sample(9)'''

'''# random test 2 sorted a ascending
# path_p:
p_xs = np.random.random_sample(11)
p_ys = np.random.random_sample(11)
p_xs.sort()
p_ys.sort()
# path_q:
q_xs = np.random.random_sample(5)
q_ys = np.random.random_sample(5)'''

'''# random test 2 sorted a and b ascending
# path_p:
p_xs = np.random.random_sample(11)
p_ys = np.random.random_sample(11)
p_xs.sort()
p_ys.sort()
# path_q:
q_xs = np.random.random_sample(11)
q_ys = np.random.random_sample(11)
q_xs.sort()
q_ys.sort()'''

'''# random test save 1
# path_p:
p_xs = [0.4772922092138784, 0.01530327285512556, 0.18623288333654608, 0.9571042784382312, 0.7531384890787559, 0.5586263445891971, 0.4087144781918828, 0.742737156990017]
p_ys = [0.9717921811331545, 0.2352575436671138, 0.5101121131543056, 0.04746020008338181, 0.45603230757057844, 0.44729409272770604, 0.4203528406682898, 0.27839889123856565]
# path_q:
q_xs = [0.46355916536642316, 0.7848048532306606, 0.8984278039692681, 0.07598172633778322, 0.2466328327882239, 0.04643325875918691, 0.14385640913682562]
q_ys = [0.2252646424367046, 0.3992611940661893, 0.1551268264371849, 0.7049183480244169, 0.7245003824862535, 0.5009416536371276, 0.2867299714276945]'''

'''# random test save 2
# path_p:
p_xs = [0.47624109293124584, 0.2837557133906694, 0.7135247957391838, 0.17614720063872424, 0.9887899709217016, 0.09772925938927324, 0.8221920560820938, 0.2368133627598178]
p_ys = [0.6917492186536913, 0.6919972280719934, 0.025707124294570338, 0.896297605839744, 0.6505998191096224, 0.5765972821537704, 0.0035238865886881854, 0.6612671752727163]
# path_q:
q_xs = [0.41887808669385074, 0.8726570820763088, 0.6650049068733719, 0.653036889162294, 0.9374457599034023, 0.8068954043422238, 0.08663305551417921]
q_ys = [0.90104618488335, 0.9961534911686879, 0.7135675484774522, 0.6274281145154397, 0.6881091124705155, 0.7559564304821617, 0.26997753895074517]'''

'''# random test save 3
# path_p:
p_xs = [0.5844585287806708, 0.946788390392606, 0.6802926446523537, 0.8527705036922676, 0.8904700476354968, 0.45888211930772094, 0.11322004405514141, 0.5373560847553378]
p_ys = [0.3036266614478819, 0.07711430842039535, 0.4951652229311967, 0.9973536237366387, 0.03730097365318219, 0.5766324346688314, 0.03316200643803002, 0.7341082819453788]
# path_q:
q_xs = [0.2714870230657075, 0.7018871765012489, 0.20877970643167865, 0.34637144332902836, 0.36469879968808705, 0.8214024752368971, 0.43857466005390366]
q_ys = [0.2990839304670133, 0.9425803453259921, 0.1871239773564659, 0.5761015296669191, 0.25768162124126626, 0.22270074761139713, 0.39085095825645777]'''

# random test save 4 - test new type of events here !!!
# path_p:
p_xs = [0.8551000828942237, 0.16531452408492242, 0.564186564823705, 0.15137331025640277, 0.9196383064638842, 0.10975251213077752]
p_ys = [0.4788005011702601, 0.13212424890044872, 0.819807938966997, 0.5166440244970916, 0.8808439598407818, 0.7282421177061701]
# path_q:
q_xs = [0.25420679369343846, 0.3463219160682829, 0.2672613971484459, 0.020871754464840686, 0.030891347542877257, 0.8818101325621397]
q_ys = [0.32744502894415317, 0.39350844910524163, 0.4760166754222921, 0.7233726864949581, 0.3406476659866665, 0.5683792989898663]

'''# random test save 5
# path_p:
p_xs = [0.6481483062554458, 0.8578194992242638, 0.8943701961742954, 0.7768759901631347, 0.6973444376862212, 0.13697759422157096]
p_ys = [0.6509430436341074, 0.34426139327156846, 0.5066826647727425, 0.5951414245061053, 0.5863777694708836, 0.6346058661988038]
# path_q:
q_xs = [0.7851525078998433, 0.2524356657434448, 0.5092419109623538, 0.4384261256283123, 0.023352609263836532, 0.3221929659258621]
q_ys = [0.6246897947866462, 0.26849238647593787, 0.21940702643448529, 0.2958678582265848, 0.4665278809112351, 0.9566606918792832]'''

'''# random test save 6
# path_p:
p_xs = [0.00599793410059446, 0.8074096648485326, 0.30171961830583205, 0.6594826159111393, 0.15829708187549696, 0.9726756942236705]
p_ys = [0.6218917474138568, 0.7125551605992374, 0.027672428048585607, 0.20596081914714048, 0.8222082488123174, 0.49302729509447085]
# path_q:
q_xs = [0.21804770554647657, 0.6437412229747035, 0.692738582481725, 0.08090341390223466, 0.8221408554840375, 0.8108620330253187]
q_ys = [0.9970005478897139, 0.6063842539417315, 0.5235865820520056, 0.18795322768660194, 0.14416811574359778, 0.06089608029252358]'''

'''# random test save 7
# path_p:
p_xs = [0.8387003590242282, 0.986386796506964, 0.9890608416994826, 0.2203996716415717, 0.29551696111733516, 0.90554071737279]
p_ys = [0.4445263246241745, 0.025871041183063825, 0.8866127110341312, 0.1596431786118372, 0.4425636036188383, 0.2042474799274986]
# path_q:
q_xs = [0.9758684134278881, 0.9718634289622968, 0.34757983101296774, 0.7876781521302637, 0.655706807877401, 0.8544511150947007]
q_ys = [0.7019170394710076, 0.23212285730847304, 0.6191001539175823, 0.8836449724340486, 0.037519589793276964, 0.2809519892742768]'''

'''# random test save 8
# path_p:
p_xs = [0.8591176539658719, 0.2920315423434664, 0.3992420367080499, 0.9145207405211361]
p_ys = [0.09890376541991241, 0.418883273171729, 0.6373631794991932, 0.15196780995368941]
# path_q:
q_xs = [0.5081380946296423, 0.6977185619019116, 0.47250825455691914]
q_ys = [0.7425394847018348, 0.2869105596549423, 0.3057053683087]'''

'''# random test save 9
# path_p:
p_xs = [0.47479273660946175, 0.16025346335529966, 0.43177138699884376, 0.9246314175816142, 0.10685165177261136,
        0.5613190486246639, 0.6920659125314877, 0.31388162717608603, 0.3651080432852446, 0.560470177868245,
        0.42823589028758346]
p_ys = [0.2742164200614492, 0.2458947011944772, 0.9910618920907522, 0.3100895358575001, 0.9858334145480506,
        0.28613757536029016, 0.7641822382983995, 0.28728153889499153, 0.9720500318213272, 0.9776421399731092,
        0.09053043938817351]
# path_q:
q_xs = [0.7653055796172689, 0.4468612051945503, 0.571562619968749, 0.46795107875269515, 0.5021143141841385,
        0.6682351819907811, 0.1834579257977117, 0.4238562175227998, 0.7009169600669711, 0.41305451875652033,
        0.5625537296401396]
q_ys = [0.08481584394654784, 0.5724664432419776, 0.7997432398252383, 0.3429615394238863, 0.46464983618715683,
        0.7461086327330805, 0.730183537249749, 0.009817595260297907, 0.6859175874417596, 0.363534187109228,
        0.19514919115926255]'''

'''# random test save 10
# path_p:
p_xs = [0.24490404366767138, 0.3171915111007536, 0.8832839237653799, 0.6555732216110922, 0.4102314395806076,
        0.21845253282665378, 0.15904657789456367, 0.2543918767181734, 0.48285570373486586, 0.24484877437328423]
p_ys = [0.7047553990190104, 0.7595345838653745, 0.3335867994909364, 0.956042535796096, 0.5755722241115714,
        0.6030177436366394, 0.5033451728548165, 0.7501239900420063, 0.41659258195972093, 0.16163362161416406]
# path_q:
q_xs = [0.17739829528634143, 0.06645717732865575, 0.42723225303552537,
        0.8075409074894082, 0.28228201053404667, 0.44492009499740137, 0.4216583490250866, 0.29645257207041564,
        0.4489301072379075]
q_ys = [0.9603745202172611, 0.725766749354519, 0.9628474677549848,
        0.60601655403887, 0.6453080452805486, 0.9233595889692813, 0.18263074027634818, 0.5986203378966767,
        0.09951086960058875]'''

'''# random test save 11
# path_p:
p_xs = [0.70229914184596653, 0.27150107271625834, 0.079395488598568154, 0.80999398366251196, 0.8376669647878412, 0.52537971263256356]
p_ys = [0.99542422754001725, 0.11021642501006046, 0.5500968467764924, 0.92402207987125762, 0.31213473174110917, 0.53900978964127899]
# path_q:
q_xs = [0.55827626447493917, 0.61265759937883912, 0.28980034706214342, 0.68097859139958572, 0.99122776678415336]
q_ys = [0.38847802876358506, 0.81362649295292611, 0.26995522166436126, 0.11922423042385721, 0.2449203374077854]'''

'''# random test save 12
# path_p:
p_xs = [0.07368875462428004, 0.15406516689603833, 0.30746312603161996, 0.41008858867809894, 0.63778232305882909, 0.69957480278830886, 0.72509868622526696, 0.91178684118443942, 0.93263929759183473, 0.96328219826322281, 0.97647714819166254]
p_ys = [0.024266148702493551, 0.02857685272755961, 0.1027948227133294, 0.21448450596528379, 0.22794146656856151, 0.34540522626400028, 0.38534936163417643, 0.49969623246955752, 0.54321916265640713, 0.85199368379579721, 0.99332148295189082]
# path_q:
q_xs = [0.0872948264706699, 0.8521084821476651, 0.55117110474391939, 0.58431552538403708, 0.70777628616953303]
q_ys = [0.15273408435141427, 0.8589243042828173, 0.25335819049509811, 0.037473744227036532, 0.89533632499926969]'''

'''# traversal test 1
# path_p:
p_xs = [0, 1, 1, 2, 2, 3, 3]
p_ys = [0, 0, 1, 1, 2, 2, 3]
# path_q:
q_xs = [0, 0, 1, 1, 2, 2, 3]
q_ys = [0, 1, 1, 2, 2, 3, 3]'''

'''# traversal test 2
# path_p:
p_xs = [0, 1, 1, 3, 3]
p_ys = [0, 0, 2, 2, 3]
# path_q:
q_xs = [0, 0, 2, 2, 3]
q_ys = [0, 1, 1, 3, 3]'''

'''# traversal test 3
# path_p:
p_xs = [0, 1, 3, 3]
p_ys = [0, 0, 2, 3]
# path_q:
q_xs = [0, 0, 2, 2, 3]
q_ys = [0, 1, 1, 3, 3]'''

'''# traversal test 4 *
# path_p:
p_xs = [0, 3, 1, 3]
p_ys = [0, 0, 2, 3]
# path_q:
q_xs = [0, 0, 2, 2, 3]
q_ys = [0, 1, 1, 3, 3]'''

'''# traversal test 5
# path_p:
p_xs = [0, 1, 3, 3]
p_ys = [0, 0, 2, 3]
# path_q:
q_xs = [0, 2, 2, 1, 3]
q_ys = [0, 2, 1, 1, 3]'''

'''# traversal test 6 *
# path_p:
p_xs = [0, 1.5, 0.5, 3.5, 3]
p_ys = [0, 0, 2, 2, 3]
# path_q:
q_xs = [0, 0, 2, 2, 3]
q_ys = [0, 1.5, 0.5, 3.5, 3]'''

'''# traversal test 7 *
# path_p:
p_xs = [-0.5, 1, 1.5, 1, 0, -1, -1, 1, 0.5, 1.5]
p_ys = [0.5, 0.5, 1, 2, 1, 2, 3, 3, 4.5, 4.5]
# path_q:
q_xs = [0, -1, 0.5, 0, -1, -1, 0, 0, 1.5, 1.5]
q_ys = [0, 1, 1, 2, 1, 2, 2, 4, 3.5, 4.5]'''

'''# traversal test 8 *
geg_epsilon_squared = 2
x1 = 0.5 * geg_epsilon_squared
x2 = math.sqrt(0.5 * (geg_epsilon_squared - 2))
# path_p:
p_xs = [0, x1, x1 - 1, x1 + 1] + [(x1+1) + x2, (x1+1) + x2]
p_ys = [0, 0, x1 + 1, x1 + 1] + [(x1+1) + -1, (x1+1) + x2]
# path_q:
q_xs = [0, 0, x1 + 1, x1 + 1] + [(x1+1) + -1, (x1+1) + x2]
q_ys = [0, x1, x1 - 1, x1 + 1] + [(x1+1) + x2, (x1+1) + x2]'''

# traversal test 9 a **!!!
# path_p:
p_xs = [0, 4, 2, 6]
p_ys = [0, 0, 6, 6]
# path_q:
q_xs = p_ys
q_ys = p_xs

'''# traversal test 9 b **
# path_p:
p_xs = [0, 3, 3, 6]
p_ys = [0, -1, 7, 6]
# path_q:
q_xs = p_ys
q_ys = p_xs'''

'''# traversal test 9 comb **
# path_p:
x = math.sqrt(8)
p_xs = [-x, 0, 4, 2, 6] + [10, 8, 12] + [16, 16] + [19, 19, 22] + [22 + x]
p_ys = [0, 0, 0, 6, 6] + [6, 12, 12] + [12, 16] + [15, 23, 22] + [22]
# path_q:
q_xs = p_ys
q_ys = p_xs'''

# traversal test 9 a complexified 1 !!!
# path_p:
p_xs = [0, 4, 3, 4, 2, 3, 2, 6]
p_ys = [0, 0, 2, 3, 3, 4, 6, 6]
# path_q:
q_xs = [0, 0, 2, 3, 3, 4, 6, 6]
q_ys = [0, 4, 3, 2, 4, 3, 2, 6]

'''# traversal test 9 a complexified 2
# path_p:
p_xs = [0, 4, 3, 4, 2, 3, 2, 6]
p_ys = [0, 0, 2, 3, 3, 4, 6, 6]
# path_q:
q_xs = p_ys
q_ys = p_xs'''

'''# old type of critical event test 1
# path_p:
p_xs = [1, 1, 1, 1]
p_ys = [1, 2, 0, 1]
# path_q:
q_xs = [0, 2]
q_ys = [1, 1]'''

print("# path_p:")
print("p_xs = " + str(list(p_xs)))
print("p_ys = " + str(list(p_ys)))
print("# path_q:")
print("q_xs = " + str(list(q_xs)))
print("q_ys = " + str(list(q_ys)))

path_p = xy_to_vectors(p_xs, p_ys)
path_q = xy_to_vectors(q_xs, q_ys)

input1 = CellMatrix(path_p, path_q, traverse=1)
print(input1)

#sample1 = input1.sample_l(-1, n_p)
#sample1 = input1.sample([4, 2.23, 1], n_p)

traversal = input1.traversals[0]
critical_epsilons = traversal.epsilons.copy()
critical_epsilons.sort()
sample1 = input1.sample(critical_epsilons[-5:], n_p)

sample_to_matplotlib(sample1, plot_cross_sections=True, plot_critical_traversals=True, plot_traversals=True, plot_3d=True)
