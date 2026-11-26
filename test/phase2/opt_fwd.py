#!/usr/bin/env python3
# Selección hacia delante de variables del BDT (26-nov-2026): parte de una lista base y, en cada paso, entrena un BDT
# por cada variable candidata añadida (opt_fwd.C, en paralelo) y se queda con la que más baja la figura de mérito
# score = error_B/error_B(base) + error_A/error_A(base). Imprime la curva error frente a número de variables.
#
#   python3 -I test/phase2/opt_fwd.py <pop> <dir> <pasos> <procesos>
#
# Necesita el entorno de highland (root en el PATH) y los árboles <dir>/slim_<pop>.root de opt_fwd_slim.C.
import subprocess, sys, os, concurrent.futures as cf

BASE = {0: "as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop".split(","),
        1: "as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop,axis,clear".split(",")}
POOL_COMMON = ("nmich5,nmichel_own,dau_maxmichel,dau_nearstop,ar_n10,ar_charge10,as_n5,as_n10,as_mean10,as_charge10,"
               "dedx_start,intpeak,isshower,cnn_trk,cnn_em,ndau_trk,ndau_shw,diry_stop,faces,gapmin").split(",")
POOL = {0: POOL_COMMON + "vtx,vc_cover,isprimary,dstop_v,dstart_v".split(","),
        1: POOL_COMMON + "slice_min,cosmiccrit".split(",")}

def train(pop, d, varlist, tag):
    slim = os.path.join(d, "slim_%d.root" % pop)
    cmd = ["root", "-l", "-b", "-q", 'test/phase2/opt_fwd.C+("%s", %d, "%s", "%s", "%s")' % (slim, pop, ",".join(varlist), d, tag)]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("RESULTADO"):
            eB = float(line.split("error B ")[1].split(" %")[0]); eA = float(line.split("error A ")[1].split(" %")[0])
            return eB, eA
    return None

def main():
    pop, d, steps, nproc = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    cur = list(BASE[pop]); pool = [v for v in POOL[pop] if v not in cur]
    b0 = train(pop, d, cur, "p%d_base" % pop)
    print("base (%d variables): error B %.2f %%, error A %.2f %%" % (len(cur), b0[0], b0[1]), flush=True)
    for step in range(1, steps + 1):
        res = {}
        with cf.ThreadPoolExecutor(nproc) as ex:
            fut = {ex.submit(train, pop, d, cur + [v], "p%d_s%d_%s" % (pop, step, v)): v for v in pool}
            for fu in cf.as_completed(fut):
                r = fu.result()
                if r: res[fut[fu]] = r
        ranked = sorted(res.items(), key=lambda kv: kv[1][0]/b0[0] + kv[1][1]/b0[1])
        print("\npaso %d: candidatas (error B, error A), de mejor a peor:" % step)
        for v, (eB, eA) in ranked:
            print("   + %-14s %5.2f %%  %5.2f %%" % (v, eB, eA))
        best = ranked[0][0]
        cur.append(best); pool.remove(best)
        print("  -> entra %s: %d variables, error B %.2f %%, error A %.2f %%" % (best, len(cur), res[best][0], res[best][1]), flush=True)
    print("\nlista final:", ",".join(cur))

main()
