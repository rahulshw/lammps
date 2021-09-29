// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Soft-core version:  Rahul Shaw
------------------------------------------------------------------------- */

#include "pair_coul_soft2.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

PairCoulCutSoft2::PairCoulCutSoft2(class LAMMPS *lmp) :Pair(lmp) {
    writedata = 0;
}


PairCoulCutSoft2::~PairCoulCutSoft2() {
    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);

        memory->destroy(cut);
        memory->destroy(scale);
    }
}

void PairCoulCutSoft2::compute(int eflag, int vflag) {
    int i,j,ii,jj,inum,jnum,itype,jtype;
    double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
    double rsq,r2inv,rinv,forcecoul,factor_coul;
    int *ilist,*jlist,*numneigh,**firstneigh;

    ecoul = 0.0;
    ev_init(eflag,vflag);

    double **x = atom->x;
    double **f = atom->f;
    double *q = atom->q;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    double *special_coul = force->special_coul;
    int newton_pair = force->newton_pair;
    double qqrd2e = force->qqrd2e;
    double rsq_by_alphasq, exp_factor;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // loop over neighbors of my atoms

    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        qtmp = q[i];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            factor_coul = special_coul[sbmask(j)];
            j &= NEIGHMASK;

            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            jtype = type[j];

            if (rsq < cutsq[itype][jtype]) {
                r2inv = 1.0/rsq;
                rinv = sqrt(r2inv);


                forcecoul = qqrd2e * scale[itype][jtype] * qtmp*q[j]*rinv;

                if (alpha[itype][jtype] > 0.0){
                    rsq_by_alphasq = rsq/(alpha[itype][jtype]*alpha[itype][jtype]);
                    exp_factor = exp(-1*rsq_by_alphasq);
                    forcecoul -= qqrd2e * scale[itype][jtype] * qtmp * q[j]
                                 * (1+2*rsq_by_alphasq) * exp_factor * rinv;
                }

                fpair = factor_coul*forcecoul * r2inv;

                f[i][0] += delx*fpair;
                f[i][1] += dely*fpair;
                f[i][2] += delz*fpair;
                if (newton_pair || j < nlocal) {
                    f[j][0] -= delx*fpair;
                    f[j][1] -= dely*fpair;
                    f[j][2] -= delz*fpair;
                }

                if (eflag) {
                    ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp * q[j] * rinv;
                    if (alpha[itype][jtype] > 0.0) {
                        ecoul -= exp_factor;
                    }
                }

                if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                     0.0,ecoul,fpair,delx,dely,delz);
            }
        }
    }

    if (vflag_fdotr) virial_fdotr_compute();
}

void PairCoulCutSoft2::allocate() {
    allocated = 1;
    int n = atom->ntypes;

    memory->create(setflag, n+1, n+1, "pair:setflag");
    memory->create(scale, n+1, n+1, "pair:scale");

    // this could have been optimized for less memory space
    // but the number of types of atoms is very less,
    // so it will not be much helpful
    for (int i=1; i<=n; i++){
        for(int j=i; j<=n; j++){
            setflag[i][j] = 0;
            scale[i][j] = 1.0;
        }
    }

    memory->create(cutsq, n+1, n+1, "pair:cutsq");
    memory->create(cut, n+1, n+1, "pair:cut");
    memory->create(alpha, n+1, n+1, "pair:alpha");
}


void PairCoulCutSoft2::settings(int narg, char **arg) {
    if (narg != 2) error->all(FLERR, "Illegal pair_style command. Expecting alpha and cutoff as arguments");


    alpha_global = utils::numeric(FLERR, arg[0], false, lmp);
    cut_global = utils::numeric(FLERR, arg[1], false, lmp);

    if (allocated) {
        int n = atom->ntypes;
        for (int i = 1; i <= n; i++) {
            for (int j = i; j<=n; j++){
                if (setflag[i][j]) {
                    cut[i][j] = cut_global;
                    alpha[i][j] = alpha_global;
                }
            }
        }
    }
}


void PairCoulCutSoft2::coeff(int narg, char **arg) {
    if (narg < 2 || narg > 4)
        error->all(FLERR, "Incorrect args for coul/cut/soft2 pair coefficients");
    if (!allocated) allocate();

    int ilo, ihi, jlo, jhi;
    int n = atom->ntypes;

    utils::bounds(FLERR, arg[0], 1, n, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, n, jlo, jhi, error);

    double alpha_one = alpha_global;
    if (narg>=3) alpha_one = utils::numeric(FLERR, arg[2], false, lmp);

    double cut_one = cut_global;
    if (narg==4) cut_one = utils::numeric(FLERR, arg[3], false, lmp);

    int count = 0;
    for (int i = ilo; i<= ihi; i++){
        for (int j = MAX(jlo, i); j<=jhi; j++){
            cut[i][j] = cut_one;
            alpha[i][j] = alpha_one;
            alpha[j][i] = alpha_one;

            scale[i][j] = 1.0;
            setflag[i][j] = 1;
            count++;
        }
    }

    if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


void PairCoulCutSoft2::init_style() {
    if (!atom->q_flag)
        error->all(FLERR,"Pair style coul/cut/soft2 requires atom attribute q");

    neighbor->request(this, instance_me);
}


double PairCoulCutSoft2::init_one(int i, int j) {
    if (setflag[i][j] == 0){
        cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
        alpha[i][j] = alpha[j][i];
        scale[i][j] = 1.0;
    }

    scale[j][i] = scale[i][j];

    return cut[i][j];
}


double PairCoulCutSoft2::single(int i, int j, int itype, int jtype,
                                double rsq, double factor_coul, double factor_lj,
                                double &fforce) {
    double r2inv, rinv, forcecoul, phicoul;

    r2inv = 1/rsq;
    rinv = sqrt(r2inv);

    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv;


    fforce = factor_coul * forcecoul * r2inv;

    phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*rinv;

    return factor_coul*phicoul;
}


void PairCoulCutSoft2::write_data(FILE *fp) {
    for (int i = 1; i <= atom->ntypes; i++)
        fprintf(fp,"%d\n",i);
}

void PairCoulCutSoft2::write_data_all(FILE *fp) {
    for (int i = 1; i <= atom->ntypes; i++)
        for (int j = i; j <= atom->ntypes; j++)
            fprintf(fp,"%d %d %g\n",i,j,cut[i][j]);
}

void PairCoulCutSoft2::write_restart(FILE *) {}

void PairCoulCutSoft2::read_restart(FILE *) {}

void PairCoulCutSoft2::write_restart_settings(FILE *) {}

void PairCoulCutSoft2::read_restart_settings(FILE *) {}

void *PairCoulCutSoft2::extract(const char *, int &) { return nullptr; }
