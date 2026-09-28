"""Numerical finite inverse seed; interval code must certify it separately."""

import argparse
import time
import numpy as np
from disc_rank2 import Disc


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--conformal',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    data=np.load(a.conformal)
    m,step,M,S=map(int,(data['m'],data['step'],data['M'],data['S']))
    scale=float(data['boundary_scale']) if 'boundary_scale' in data else 1.0
    D=Disc(m,M,S,step=step,boundary_scale=scale)
    start=time.monotonic()
    J=D.jacobian(data['z'])
    print('matrix',J.shape,'seconds',time.monotonic()-start,flush=True)
    A=np.linalg.inv(J)
    print('inverse seconds',time.monotonic()-start,
          'floating defect',np.linalg.norm(np.eye(len(J))-A@J,ord=np.inf),flush=True)
    np.save(a.output,A)


if __name__=='__main__':
    main()
