from lv_init import *
import json
import glob
from findstr import *

#------------------------------------------------- INITIALIZATION ------------------------------------------------------
savefolder = 'mto_opt_2'
savepath = os.path.dirname(path) + '/' + savefolder + '/'
fit_folder = 'fitness'
target_folder = 'targets'
cons_folder = 'constraints'
out_folder = 'outputs'
target_path = savepath + target_folder + '/'
fit_path = savepath + fit_folder + '/'
cons_path = savepath + cons_folder + '/'
out_path = savepath + out_folder + '/'
ngen_start = 1

pop_size = 160
ngen = 10000
k0 = 1
k1 = 1
cross_prob = 0.7
mut_prob = 0.5
# Initialization for Optimizer
f_sum = np.array([0, 0, 0])
v_stage = np.array([0, 0, 0])
mdot = np.array([0, 0, 0])
sf = strf
# Isp = Sp_imp
bnds = [[0, 0]]*11
bnds[0] = np.array([8e3, 11e3])
bnds[1] = np.radians([-2, 0])
bnds[2] = np.radians([-2, 2])
bnds[3] = np.radians([-2, 2])
bnds[4] = np.array([2e3*g0, 40e3*g0])
bnds[5] = np.array([0.5e3*g0, 10e3*g0])
bnds[6] = np.array([0.01, 0.99])
bnds[7] = np.array([0.01, 0.99])
bnds[8] = np.array([0.01, 0.99])
bnds[9] = np.array([0.01, 0.99])
bnds[10] = np.array([0.01, 0.99])
segp = np.array([1, 5, 5])
segt = np.array([1, 5, 5])
sl = []
sl.extend(segp.tolist()), sl.extend(segt[1:].tolist())
seg_len = [1]+sl*2
# norm_list = [bnds[0][0], -bnds[1][0], bnds[2][1], bnds[3][1], bnds[4][1], ]

if not os.path.exists(savepath):
    os.makedirs(savepath)
if not os.path.exists(target_path):
    os.makedirs(target_path)
if not os.path.exists(fit_path):
    os.makedirs(fit_path)
if not os.path.exists(cons_path):
    os.makedirs(cons_path)
if not os.path.exists(out_path):
    os.makedirs(out_path)

H_req = 500e3
V_req = 1000 * np.sqrt(398600 / (Rs / 1000 + 500))
Theta_req = 0

feas_alt = 5e3 # 10 kilometres
feas_fpa = 1*np.pi/180 # 2 degrees
feas_acc = 4*g0 # 4g acceleration
feas_vel = 8 # 50 m/s

if ngen_start > 1:
    target_list = glob.glob(target_path + '*')
    sorted_list = sorted(target_list, key=os.path.getmtime)
    tar_filename = sorted_list[-3]
    # tar_filename = target_path + '4015_gen_targets'

    n_start = findOccurrences(sorted_list[-1], '/')[-1]
    n_end = findOccurrences(sorted_list[-1], '_')[-2]
    ngen_start = int(sorted_list[-3][n_start + 1:n_end])+1
    # ngen_start = 4016

    with open(tar_filename, 'r') as x:
        inputs = json.load(x)

    print('\nLoading population from Target file:', tar_filename, '\n')

    targets = []
    for tar in inputs:
        target = []
        for inp in tar[0]:
            target.append(np.array(inp))
        targets.append([target])

    pop_size = len(targets)  # population size

else:

    print('\nTarget file not readable. Creating random population.\n')
    targets = []
    for i in range(pop_size):
        np.random.seed(i); tar0 = np.random.uniform(bnds[0][0], bnds[0][1], 1) # Ideal Velocity
        np.random.seed(i); tar1 = np.random.uniform(bnds[1][0], bnds[1][1], segp[0]) # Kickrate
        np.random.seed(i); tar2 = np.random.uniform(bnds[2][0], bnds[2][1], segp[1]) # Pitchrates, 2nd stage
        np.random.seed(i); tar3 = np.random.uniform(bnds[3][0], bnds[3][1], segp[2]) # Pitchrates, 3rd stage
        np.random.seed(i); tar4 = np.random.uniform(bnds[4][0], bnds[4][1], segt[1])  # Thrusts, 2nd stage
        np.random.seed(i); tar5 = np.random.uniform(bnds[5][0], bnds[5][1], segt[2])  # Thrusts, 3rd stage
        np.random.seed(i); tar6 = np.random.rand(segp[0]) # Kickrate interval, 1st stage
        np.random.seed(i); tar7 = np.random.rand(segp[1]) # Pitchrate intervals, 2nd stage
        tar7 = tar7 / np.sum(tar7)
        np.random.seed(i); tar8 = np.random.rand(segp[2]) # Pitchrate intervals, 3rd stage
        tar8 = tar8 / np.sum(tar8)
        np.random.seed(i); tar9 = np.random.rand(segt[1]) # Thrust intervals, 2nd stage
        tar9 = tar9 / np.sum(tar9)
        np.random.seed(i); tar10 = np.random.rand(segt[2]) # Thrust intervals, 3rd stage
        tar10 = tar10 / np.sum(tar10)

        target = [[tar0, tar1, tar2, tar3, tar4, tar5, tar6, tar7, tar8, tar9, tar10]]
        targets.append(target)

input_sav = {
	     	'Np': pop_size,
		'F': mut_prob,
		'Cr': cross_prob,
        'Isp': Sp_imp.tolist(),
        'sf': strf.tolist(),
        'Thrust': Thr.tolist(),
        'Thrust Segments' : segt.tolist(),
        'Pitchrate Segments' : segp.tolist(),
        'Left out prop': m_lo.tolist(),
        'Altitude Constraint': feas_alt,
        'Injection Angle Constraint' : feas_fpa,
        'Acceleration Constraint': feas_acc,
		'folder name': savefolder,
        'Comments' : 'Lift off Mass Optimization'
	    }

y = json.dumps(input_sav)
f = open(savepath + 'input_sav','w')
f.write(y)
f.close()

fit_save1 = []
fit_save2 = []
fit_save3 = []
tar_norm_save = []

gen_no = ngen_start

def list_to_seg(lists):
    seg_out = []
    for lis in lists:
        si = 0
        seg = []
        for ite in seg_len:
            seg.append(lis[si:si+ite])
            si+=ite
        seg_out.append([seg])
    return seg_out


def seg_to_list(segs):
    lists_out = []
    for segx in segs:
        list_out = []
        for sx in segx[0]:
            list_out.extend(sx)
        lists_out.append(list_out)
    return lists_out
