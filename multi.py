from pipe import Pipeline
from multiprocessing import Process

if __name__ == "__main__":
    p1 = Pipeline()
    p2 = Pipeline(circular=True)
    p3 = Pipeline(clocklike=True)
    p4 = Pipeline(circular=True, clocklike=True)

    P  = [p1, p2, p3, p4]

    N = 10
    min_leafes = 6
    max_leafes = 9
    first_leafes = [1,2,3,4]
    B = [0,1,2,3]
    leaves_unknown = True
    subset = 4
    choose_random = True
    small_spike = True
    seed_szenario = 0
    seed_pipeline = 0
    seed_permutations = 0

    L:list[Process] = []
    for p in P:
        l1 = Process(target=p.wp1, args=(N, min_leafes, max_leafes, first_leafes, 
            seed_szenario, "", ))
        l2 = Process(target=p.wp2, args=(N, min_leafes, max_leafes, first_leafes,
            choose_random, seed_szenario, seed_pipeline, "", ))
        l3 = Process(target=p.wp3, args=(N, min_leafes, max_leafes, first_leafes, 
            B, leaves_unknown, subset, False, seed_szenario, seed_pipeline, 
            seed_permutations, "", ))
        l4 = Process(target=p.wp4, args=(N, min_leafes, max_leafes, first_leafes,
            small_spike, True, seed_szenario, seed_pipeline, ""))
            
        L += [l1, l2, l3, l4]

    for l in L:
        l.start()
    
    for l in L:
        l.join()