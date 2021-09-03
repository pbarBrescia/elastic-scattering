from hyperopt import fmin, tpe
from hyperopt import hp

from objective import obj

space = hp.choice()

best = fmin(obj, space, algo=tpe.suggest, max_evals=100)