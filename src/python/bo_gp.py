import numpy as np
import time
import os
import re

import matplotlib.pyplot as plt

from skopt import Optimizer
from skopt import Space
from skopt import dump, load

from skopt.plots import plot_objective, plot_convergence, plot_evaluations

from src.python.objective import obj

np.random.seed(1234)


class BoRichGp:
    def __init__(self, obj, space, dir='.', n_call=10, id="test"):
        """
        Class for applying bayesian optimization.

        :param obj:
        :param space:
        :param dir:
        :param n_call:
        :param id:

        """
        self.obj = obj
        self.dir = dir
        self.n_call = n_call
        self.space = space
        self.id = id

        self.opt = Optimizer(space,
                "GP",
                n_initial_points=3,
                acq_optimizer="sampling",
                initial_point_generator='lhs')

    def optimize(self, keep=True):
        old, file, number = self._chek_last_checkpoint()
        if old:
            self.load_last_checkpoint(file)

        start = time.time()
        for i in range(int(self.n_call)):
            x = self.opt.ask()
            y = self.obj(x)
            self.opt.tell(x, y)
            if keep:
                dump(self.opt.get_result(), self.dir + '/' + f'checkpoint_{self.id}_{number}')
        stop = time.time()
        print(f'time spent for {self.n_call} calls of {self.id}: ', stop-start, ' seconds')

    def plot_conv(self):
        _ = plot_convergence(self.opt.get_result())
        plt.savefig(self.dir + '/' + f'conv_{self.id}.png')
        plt.show()

    def plot_depend(self):
        _ = plot_objective(self.opt.get_result())
        plt.savefig(self.dir + '/' + f'obj_{self.id}.png')
        plt.show()

    def plot_eval(self):
        _ = plot_evaluations(self.opt.get_result())
        plt.savefig(self.dir + '/' + f'eval_{self.id}.png')
        plt.show()

    def _chek_last_checkpoint(self):
        file_number = 0
        check = False
        file_name = ""
        regex = re.compile(r'\d+')

        for file in os.listdir(self.dir):
            if ('checkpoint_' + self.id) in file:
                cur_file_number = int(regex.findall(file)[0])
                if cur_file_number >= file_number:
                    file_number = cur_file_number
                    file_name = file
                check = True

        return check, file_name, file_number + 1

    def load_last_checkpoint(self, file_name):
        res = load(self.dir + '/' + file_name)
        x0 = res.x_iters
        y0 = res.func_vals
        for point in range(len(x0)):
            self.opt.tell(x0[point], y0[point])

        return res


if __name__ == "__main__":

    SPACE = Space.from_yaml('space.yaml')
    opt = BoRichGp(obj=obj, space=SPACE, id='namedim_test', n_call=40)

    opt.optimize()
    opt.plot_conv()

    # opt.load_last_checkpoint('checkpoint_namedim_test_3')
    # opt.plot_depend()

###############################
#BO GP WITH SCIKIT
###############################
#
# noise_level = 0.1

# # Space of parameters
# SPACE = [
#     Real(0, 500, name='W'),
#     Real(0, 500, name='V'),
#     Real(0.1, 2, name='r1'),
#     Real(0.1, 2, name='r2'),
#     Real(0.1, 1, name='a1'),
#     Real(0.1, 1, name='a2')
# ]
# import re
#
# file_number = 0
# x0 = None
# y0 = None
#
# regex = re.compile(r'\d+')
# for file in os.listdir("."):
#     if file.endswith(".pkl"):
#         cur_file_number = int(regex.findall(file)[0])
#         if cur_file_number > file_number:
#             file_number = cur_file_number
#
# if file_number > 0:
#     filename = f'checkpoint_{file_number}.pkl'
#     res = load(filename)
#     x0 = res.x_iters
#     y0 = res.func_vals
#
# checkpoint_saver = CheckpointSaver(f"./checkpoint_{file_number + 1}.pkl", compress=9) # keyword arguments will be passed to `skopt.dump`
#
# # Here the magic happens
#
# start = time.time()
#
# result = gp_minimize(obj,
#                      SPACE,
#                      n_calls=50,
#                      n_jobs=-1,
#                      x0=x0,
#                      y0=y0,
#                      noise=noise_level,
#                      acq_func="LCB",
#                      callback=[checkpoint_saver])
#
# stop = time.time()
#
# print(result)
#
#
###############################
#BO GP WITH SCIKIT
###############################

# noise_level = 0.1

# # Space of parameters
# SPACE = [
#     Real(0, 500, name='W'),
#     Real(0, 500, name='V'),
#     Real(0.1, 2, name='r1'),
#     Real(0.1, 2, name='r2'),
#     Real(0.1, 1, name='a1'),
#     Real(0.1, 1, name='a2')
# ]
# import re
#
# file_number = 0
# x0 = None
# y0 = None
#
# regex = re.compile(r'\d+')
# for file in os.listdir("."):
#     if file.endswith(".pkl"):
#         cur_file_number = int(regex.findall(file)[0])
#         if cur_file_number > file_number:
#             file_number = cur_file_number
#
# if file_number > 0:
#     filename = f'checkpoint_{file_number}.pkl'
#     res = load(filename)
#     x0 = res.x_iters
#     y0 = res.func_vals
#
# checkpoint_saver = CheckpointSaver(f"./checkpoint_{file_number + 1}.pkl", compress=9) # keyword arguments will be passed to `skopt.dump`
#
# # Here the magic happens
#
# start = time.time()
#
# result = gp_minimize(obj,
#                      SPACE,
#                      n_calls=50,
#                      n_jobs=-1,
#                      x0=x0,
#                      y0=y0,
#                      noise=noise_level,
#                      acq_func="LCB",
#                      callback=[checkpoint_saver])
#
# stop = time.time()
#
# print(result)
#
#
# print('=================================')
# print('time: ', stop-start)
#
# _ = plot_objective(result, dimensions=['W','V','r1','r2','a1','a2'])
# plt.savefig(f'obj_{file_number + 1}.png')
# plt.show()
# _ = plot_convergence(result, yscale='log')
# plt.savefig( f'conv_{file_number + 1}.png')
# plt.show()
# _ = plot_evaluations(result)
# plt.savefig(f'eval_{file_number + 1}.png')
# plt.show()

##PROVA DI DUMPING CON OPTIMIZER CLASS

# from joblib import Parallel, delayed
# from skopt import dump, load
#
#
# init_poit = 20
# n_call = 100
# n_core = 5
# file_number = 0
#
# regex = re.compile(r'\d+')
# for file in os.listdir("."):
#     if 'bunch_trial' in file:
#         cur_file_number = int(regex.findall(file)[0])
#         if cur_file_number > file_number:
#             file_number = cur_file_number
#
# if file_number > 0:
#     filename = f'bunch_trial_{file_number}.z'
#     res = load(filename)
#     x0 = res.x_iters
#     y0 = res.func_vals
#
# from skopt import Optimizer
#
# opt = Optimizer([(40., 100.),
#                  (40., 100.),
#                  (0.8, 1.5),
#                  (0.8, 1.5),
#                  (0.3, 0.8),
#                  (0.3, 0.8)],
#                 "GP",
#                 n_initial_points=init_poit,
#                 # acq_func="EI",
#                 acq_optimizer="sampling",
#                 initial_point_generator='lhs')
#
#
# start = time.time()
#
# for i in range(int(n_call/n_core)):
#     x = opt.ask(n_points=n_core)
#     y = Parallel(n_jobs=n_core)(delayed(obj)(v) for v in x)
#     opt.tell(x, y)
#     dump(opt.get_result(), f'bunch_trial_{file_number + 1 + i}.z')
#
# stop = time.time()
# # opt.run(objective, n_iter=50)
#
# result = opt.get_result()
# print(result)
#
# print('=================================')
# print('time: ', stop-start)
#
# _ = plot_objective(result, dimensions=['W','V','r1','r2','a1','a2'])
# plt.savefig(f'obj_{file_number + 1}.png')
# plt.show()
# _ = plot_convergence(result, yscale='log')
# plt.savefig( f'conv_{file_number + 1}.png')
# plt.show()
# _ = plot_evaluations(result)
# plt.savefig(f'eval_{file_number + 1}.png')
# plt.show()
#

# [99, 41, 0.9328003728491479, 1.4035506822647186, 0.4229878383297476, 0.3756859417357024]

###############################
#BO WTIH PYTORCH (using CUDA)
###############################


# class BayesianOptimization(object):
#     def __init__(self, gp_model, acq, lr=1e-1, ndim=1, n_multi_start=32, warmup=2):
#         self.gp_model = gp_model
#         self.gp = None
#
#         self.acq = acq
#         self.ndim = ndim
#         self.n_multi_start = n_multi_start
#         self.warmup = warmup
#
#         self.X = list()
#         self.y = list()
#
#     def reset(self, ):
#         self.X = list()
#         self.y = list()
#         self.gp = None
#
#     def ask(self):
#         if len(self.X) < self.warmup:
#             return np.random.uniform(0, 1, size=self.ndim)
#
#         ### optimizing aquisition function
#         suggestions, values = minimize_acq(
#             x0=np.random.uniform(0, 1, size=(self.n_multi_start, self.ndim)),
#             gp=self.gp,
#             acq=self.acq
#         )
#
#         ### returning the best guess
#         best = np.argmin(values)
#         return suggestions[best]
#
#     def tell(self, x, y):
#         self.X.append(x)
#         self.y.append(y)
#
#         X = torch.tensor(self.X, dtype=torch.float32, device=DEVICE).view(-1, self.ndim)
#         y = torch.tensor(self.y, dtype=torch.float32, device=DEVICE).view(-1, )
#
#         ### 'retraining' model
#         self.gp = self.gp_model(X, y)
#
# def minimize_acq(x0, gp, acq, lr=1e-1, n_iters=128, progress_bar=lambda x: x):
#     x = torch.tensor(x0, dtype=torch.float32, device=DEVICE, requires_grad=True)
#
#     ### not stochastic in this case
#     opt = torch.optim.SGD(lr=lr, params=[x])
#
#     for _ in progress_bar(range(n_iters)):
#         opt.zero_grad()
#         torch.sum(acq(gp, x)).backward()
#         opt.step()
#
#         with torch.no_grad():
#             x = torch.clamp(x, 0, 1)
#
#     with torch.no_grad():
#         values = acq(gp, x)
#
#     return x.detach().cpu().numpy(), values.detach().cpu().numpy()