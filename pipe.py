# -*- coding: utf-8 -*-

import os, logging, time

from erdbeermet.simulation import simulate, Scenario
from erdbeermet.recognition import recognize, is_pseudometric
from erdbeermet.tools.Tree import Tree, TreeNode

import numpy as np
from numpy.random import Generator, PCG64
import pandas as pd

from itertools import permutations


# ---------------------------
class Pipeline:
    """Create a Pipeline for all the simulations in 
    workpackages (see 'praktikumsanleitung')

    Attributes
    ----------
    circular : bool, optional
        Create only R-Matrix 'historys' with 'erdbeermet.simulation' which are
        guaranteed to be 'circular' Type R-matrices. This option will be used 
        for all simulations in a specific pipeline. The default is False.
    clocklike : bool, optional
        Create only R-Matrix 'historys' with 'erdbeermet.simulation' in which
        all distance increment are equal for all items within each iteration.
        The default is False.
    first_candidate_only : bool, optional
        Use this for the recognition algorithm in the erdbeermet package 
        'erdbeermet.recognition' as option. If True, only the first found
        candidate will be considered for the merge event. The default is
        True.

    See also
    --------
    erdbeermet.simulation
    erdbeermet.recognition
    
    """

    def __init__(self, circular:bool=False, clocklike:bool=False, first_candidate_only:bool=True):
        self.circular = circular
        self.clocklike = clocklike
        self.first_candidate_only = first_candidate_only


    # @property
    # def first_candidate_only(self):
    #     return self.__first_candidate_only

    # @first_candidate_only.setter
    # def first_candidate_only(self, first_candidate_only:bool):
    #     try:
    #         if not isinstance(first_candidate_only, bool):
    #             print("[Note] Attribute 'first_candidate_only' must be boolean.", end=" ")
    #             print(f"      Nothing changed (first_candidate_only={self.first_candidate_only}).")
    #         else:
    #             print(f"Call setter, {first_candidate_only}")
    #             self.__first_candidate_only = first_candidate_only

    #     except AttributeError:
    #         print(f"[Error] Class attribute not set correctly. (cicular = '{first_candidate_only}')")
    #         print(f"        Attribute 'first_candidate_only' is set to 'True'.")
    #         self.__first_candidate_only = True


    @property
    def circular(self):
        return self.__circular

    @circular.setter
    def circular(self, circular:bool):
        try:
            if not isinstance(circular, bool):
                print("[Note] Attribute 'circular' must be boolean.", end=" ")
                print(f"      Nothing changed (circular={self.circular}).")
            else:  
                self.__circular = circular

        except AttributeError:
            print(f"[Error] Class attribute not set correctly. (cicular = '{circular}')")
            print(f"        Attribute 'circular' is set to 'False'.")
            self.__circular = False


    @property
    def clocklike(self):
        return self.__clocklike

    @clocklike.setter
    def clocklike(self, clocklike:bool):
        try:
            if not isinstance(clocklike, bool):
                print("[Note] Attribute 'clocklike' must be boolean.", end=" ")
                print(f"      Nothing changed (clocklike={self.clocklike}).")
            else:    
                self.__clocklike = clocklike

        except AttributeError:
            print(f"[Error] Class attribute not set correctly. (clocklike = '{clocklike}')")
            print(f"        Attribute 'clocklike' is set to 'False'.")
            self.__clocklike = False


    @staticmethod
    def _find_common_steps(tree:Tree, history:list):
        """Find common R-Steps in recognition tree and history. Help function for
        '_analyse_scenario'.

        Parameter
        ---------
        tree : Tree
            Recognition Tree from function 'recognize' in 'redbeermet.recognize'.
        history : list of tuples
            List of tuples generated from Class 'Simulation' in 'erdbeermet.simulation'.
            Tuples have the format (x, y, z, u, alpha).

        Returns
        -------
        valid_ways : list
            List of nodes, which are correct recognize as R-Matrices.
        r_steps_order : int
            Number of common R-Steps in correct order in recognition tree and history.
        r_steps : int
            Number of common R-Steps recognition tree and history.
        """

        # Extract simulation steps from history ('true' triples x, y, z)
        # swap x and y if x > y (because 'recognition tree' have it)
        true_steps = []
        for i in range(len(history)):
            step = list(history[i])[:3]
            if step[0] < step[1]:
                true_steps.append([step[0], step[1], step[2]])
            else:
                true_steps.append([step[1], step[0], step[2]])
        
        valid_ways = []
        r_steps_order = 0
        r_steps = 0

        # traverse tree in preorder
        for node in tree.preorder():
            if node.n == 4 and node.valid_ways != 0:
                valid_ways.append(node)
            
            # count common steps (reconstruction vs. simulation)
            if node.R_step is not None and node.valid_ways != 0:
                if list(node.R_step)[:3] in true_steps and len(true_steps) > 3:
                    r_step = list(node.R_step)[:3]

                    if r_step == true_steps[-1]:
                        # Rstep in order == last element in list
                        r_steps_order += 1
                        r_steps += 1 
                        true_steps.pop()
                    else:
                        r_steps += 1
                        true_steps.pop()
        
        return valid_ways, r_steps_order, r_steps


    def _analyse_scenario(self, scenario:Scenario, B:list, first_leaves:list, 
        choose_random:bool, small_spike:bool, generator:Generator, info:bool):
        """Pipeline for classification process of distance matrices.

        Parameter
        ---------
        scenario : class 'Scenario'
            Generated Class 'Scenario' from the 'erdbeermet' package. Includes
            the distance matrix and history of merge and branching events.
        B: list
            Set of leaf-identifiers. If not empty no leaf of this set will we 
            chosen as z in the recognize algorithmen (see also 'erdbeermet.recognize')
        first_leaves: list
            First four leaves of the generation process.
        choose_random : bool
            From the list of correctly determined R-Maps, randomly select one 
            and evaluate whether it is the first candicate.
        small_spike : bool
            If True, use recognition algorithm by smallest spike consumption for number of 
            leaves > 5.
        generator : Generator
            Random number generator from numpy.random to make results reproducible.
        info : bool
            Print some infos from 'erdbeermet.recognize' function 'recognize'.

        Returns
        -------
        Tree: class
            The recognition tree (see also 'erdbeermet.tools.Tree')
        Dictionary:
            Result summary
        """

        D = scenario.D
        if not is_pseudometric(D):
            print("[Error] Scenario cannot be analyzed. The given matrix is not a pseudo metric. Failure!")
            return None, None

        tree = recognize(D=D, B=B, first_candidate_only=self.first_candidate_only, 
            small_spike=small_spike, print_info=info)

        # check root if recognition tree has only dead ends
        complete_dead = True if tree.root.valid_ways == 0 else False

        # Some variables
        r_metrics:list[TreeNode] = [] # list of nodes from type 'TreeNode'
        number_r_metrics = 0
        find_first_leaves = False
        select_first_leaves = False
        r_steps_order = 0
        r_steps = 0

        if not complete_dead:
            r_metrics, r_steps_order, r_steps = self._find_common_steps(tree=tree, history=scenario.history)
            number_r_metrics = len(r_metrics)

            # check if the first leaves are under the recognized pseudometrics
            for metric in r_metrics:
                if metric.V == first_leaves:
                    find_first_leaves = True
                    select_first_leaves = True
                    break
            
            # choose winner candidate at random from final r metrics
            if choose_random:
                winner_index = 0
                if number_r_metrics > 1:
                    winner_index = generator.choice(len(r_metrics)-1, size=1)

                winner_r_metric = r_metrics.pop(int(winner_index))
                for loser in r_metrics:
                    loser.valid_ways = 0
                    loser.info = "random loser"
                
                # check if we choose our first four leaves as winner
                if winner_r_metric.V != first_leaves:
                    select_first_leaves = False
        
        summary = {
            'complete_dead': complete_dead,
            'all_r_maps': number_r_metrics,
            'find_first_4leaves': find_first_leaves,
            'select_first_4leaves': select_first_leaves,
            'common_r-steps_order': r_steps_order,
            'common_r-steps': r_steps
        }

        return tree, summary


    @staticmethod
    def analyse(scenario, first_candidate_only:bool=True, B:list=[], first_leaves:list=[0,1,2,3], permu:bool=False,
        subset:int=4, choose_random:bool=False, small_spike:bool=False, seed_random:int=None, seed_permu:int=None, 
        info:bool=False):
        """ Analyse one single simulation given by 'erdbeermet.simulation' with
        pipeline method '_analyse_scenario'.

        Parameter
        ---------
        scenario : class 'Scenario'
            A scenario given by the 'erdbeermet'-Packages. See also 'erdbeermet.simulate'.
        first_candiate_only : bool, optional
            If True, only consider the first found candidate for a merge event inside the
            function 'recognize' from the 'erdbeermet'-Package. The default is True. 
            See also 'erdbeermet.recognition'
        B : list, optional
            Leaf-identifier for the recognition algorithm. During the recognition no leaf
            from this list will be chosen to reduce as z. The default is [] (every leaf can
            be chosen). See also function '_find_candidates' from 'erdbeermet.recognition'.
        first_leafes: list, optional
            The first leaves of the history simulation. In Default case they will be [1,2,3,4].
        permu : bool, optional
            If set to True algorithm will iterate over permutations (in random order) in range of 
            k-leaves (k is equal to the current number of items) as leaf identifier B of size 
            'subset' (3 or 4) until an R-Map is correctly identified (reconition tree has no dead 
            end). In this case the result .csv file will have the following additional columns:
            'permu_count' (number of iteratios needed) 'max_permu' (maximum number of possible
            permutations) and 'permutation' (permutation as leaf_identifier B which produces the
            first R-Maps). The default is False.
        subset : int, optional
            Number of core leaves from which random permutations are created and used as 
            leav-identifier. Only useful if 'permu=True'. Default is 4.
        choose_random : bool, optional
            If set to True, in the pipeline algorithm out of the list of candidates for 
            the last step, only one random candidate (with 4 leafes) is output as a positive
            R-Metric. All other will be set to 'random loser'. The default is False.
        small_spike : bool, optional
            If True, use recognition algorithm by smallest spike consumption for number of 
            leaves > 5. The default is False.
        seed_random : int, optional
            Seed to make all randomly choosen positive candidate for the same number of N simulations
            (with same seed_szenario) reproducible. Only useful with 'choose_random=True'. By 
            default it is set to None (use a random seed between 0 and max int32).
        seed_permu : int, optional
            Seed to make all randomly choosen permutation reproducible. Only useful with 
            'permu=True'. By default it is set to None (use a random seed between 0 and max int32).
        info : bool, optional
            Print more information in the recognition algorithm.

        Returns
        -------
        Tree: class
            The recognition tree (see also 'erdbeermet.tools.Tree')
        Dictionary:
            Result summary
        
        See also
        --------
        erdbeermet.recognition
        Pipeline._analyse_scenario
        """

        if not scenario or not isinstance(scenario, Scenario):
            print("[Error] No instance of the class 'Scenario' was passed. Return 'None'.")
            return None, {}
        
        if seed_random is None:
            generator_pipe = Generator(PCG64())
        else:
            generator_pipe = Generator(PCG64(seed_random))
        
        pipe = Pipeline(first_candidate_only=first_candidate_only)
        

        if permu:
            leaves = scenario.D.shape[0]
            if seed_permu is None:
                generator_permu = Generator(PCG64())
            else:
                generator_permu = Generator(PCG64(seed_random))
            
            tree, d, permu_count, max_permu, B = pipe._wp3_permu(scenario=scenario, k_leaves=leaves,
                subset=subset, first_leaves=first_leaves, choose_random=choose_random,
                small_spike=small_spike, generator_permutation=generator_permu, 
                generator_pipeline=generator_pipe)
            d["permu_count"] = permu_count
            d["max_permu"] = max_permu
            d["permutation"] = B
            return tree, d
        
        if small_spike:
            tree, d = pipe._analyse_scenario(scenario=scenario, B=B, first_leaves=first_leaves, 
                choose_random=choose_random, small_spike=small_spike, generator=generator_pipe,
                info=info)
            cycle = pipe._wp4_traverse(tree)
            d["cycle_detected"] = cycle
            return tree, d

        return pipe._analyse_scenario(scenario=scenario, B=B, first_leaves=first_leaves, 
            choose_random=choose_random, small_spike=small_spike, generator=generator_pipe,
            info=info)
    

    @staticmethod
    def _create_logger(name:str):
        """Create logger for tracking the actual simulation process of workpackages
        """

        # create logger
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)

        # create file handler and set level
        ch = logging.FileHandler(name)
        ch.setLevel(logging.DEBUG)

        # create formatter
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # add formatter to ch
        ch.setFormatter(formatter)

        # add ch to logger
        logger.addHandler(ch)
        return logger
    

    def _wp3_permu(self, scenario:Scenario, k_leaves:int, subset:int, first_leaves:list, choose_random:bool,
        small_spike:bool, generator_permutation:Generator, generator_pipeline:Generator):
        """ Helper function for the workpackage nr. 3. Iterate over all permutations
        in range of k-leaves as leaf identifier (will not be chosen) of size 'subset' 
        (3 or 4) until an R-Map is correctly identified (reconition tree has no dead end).
        """

        subset = subset
        permu = list(permutations(range(k_leaves), subset))
        max_permu = len(permu)
        generator_permutation.shuffle(permu)
        
        permu_count = 0
        for B in permu:
            B = list(B)
            tree, d = self._analyse_scenario(scenario=scenario, B=B, first_leaves=first_leaves,
                choose_random=choose_random, small_spike=small_spike, generator=generator_pipeline, info=False)
            
            permu_count += 1
            if not d["complete_dead"]:
                return tree, d, permu_count, max_permu, B
        
        return tree, d, permu_count, max_permu, B
    
    @staticmethod
    def _wp4_traverse(tree:Tree):

        node:TreeNode = None
        for node in tree.postorder():
            if node.info == "cycle detected":
                return True
        
        return False

    def _workpackages(self, simulation:str, N:int, min_leaves:int, max_leaves:int, first_leaves:list, B:list, 
        choose_random:bool, wp3_permu:bool, subset:int, small_spike:bool, seed_szenario:int, 
        seed_pipeline:int, seed_permutations:int, add_name:str):
        """Pipeline for workpackages. Parameters will be chosen by the corresponding 
        workpackages. See Pipeline.wp1, *.wp2, *.wp3.
        """

        result_dir = f"{simulation}_results"
        if add_name:
            result_dir = result_dir + add_name
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)

        opt_sce = f"result_circular-{self.circular}_clocklike-{self.clocklike}"
        # opt_rng = f"seeds-{seed_szenario}-{seed_pipeline}"
        result_name = f"{simulation}_{opt_sce}"
        
        if add_name:
            result_name = result_name + add_name
        
        file_name = result_name
        result_name += ".csv"
        target = os.path.join(result_dir, result_name)

        # Create logger
        log_name = f"simulation_{simulation}{add_name}_" \
                   f"circular-{self.circular}_clocklike-{self.clocklike}.log"
        if add_name:
            log_name = log_name + add_name
            
        if os.path.exists(log_name):
            os.remove(log_name)
        logger = self._create_logger(log_name)
        logger.info(f"Start {result_name}")

        # Set seed
        rng_szenario = Generator(PCG64(seed_szenario))
        rng_pipeline = Generator(PCG64(seed_pipeline))
        rng_permutat = Generator(PCG64(seed_permutations))

        # Parameter and matrices
        r = max_leaves - min_leaves + 1

        column_names = ['nr', 'leaves', 'circular', 'clocklike', 'dead_end', 'all_r_maps', 
            'find_first_4leaves', 'select_first_4leaves', 'r_steps', 'r_steps_order']
        
        if wp3_permu:
            column_names += ['permu_count', 'max_permu', 'permutation']
        
        if small_spike:
            column_names += ["cycle_detected"]
        
        column_names = column_names + ['runtime']

        df = pd.DataFrame(np.nan, index=list(range(r*N)), columns=column_names)

        c = 0
        l = (r*N)/10
        start_global = time.time()
        for k in range(min_leaves, max_leaves+1):
            for i in range(N):
                scenario = simulate(k, branching_prob=0.0, circular=self.circular, 
                    clocklike=self.clocklike, generator=rng_szenario)
                
                # WP 1 / 2 / part of 3
                start = time.time()
                if wp3_permu:
                    tree, d, n_permu, max_permu, b = self._wp3_permu(scenario=scenario, k_leaves=k,
                        subset=subset, first_leaves=first_leaves, choose_random=choose_random, 
                        small_spike=small_spike, generator_permutation=rng_permutat, generator_pipeline=rng_pipeline)
                   
                elif small_spike:
                    cycle = False
                    tree, d = self._analyse_scenario(scenario=scenario, B=B, first_leaves=first_leaves,
                        choose_random=choose_random, small_spike=small_spike, generator=rng_pipeline, 
                        info=False)

                else:
                    tree, d = self._analyse_scenario(scenario=scenario, B=B, first_leaves=first_leaves,
                        choose_random=choose_random, small_spike=small_spike, generator=rng_pipeline, 
                        info=False)

                end = time.time() 

                if d['complete_dead']:
                    i_ = str(i).zfill(len(str(N)))
                    if wp3_permu:
                        name_dead = f"dead_end_{file_name}_seeds" + \
                                f"-{seed_szenario}-{seed_pipeline}-{seed_permutations}" + \
                                f"_k-{k}_i-{i_}"
                    else:
                        name_dead = f"dead_end_{file_name}_seeds" + \
                                f"-{seed_szenario}-{seed_pipeline}" + \
                                f"_k-{k}_i-{i_}"

                    if small_spike:
                        cycle = self._wp4_traverse(tree)
                        name_dead = name_dead + "_cycle_" + str(cycle)

                    scenario.write_history(os.path.join(result_dir, name_dead))

                res = [
                    int(i),
                    int(k), 
                    str(self.circular), 
                    str(self.clocklike), 
                    str(d['complete_dead']),
                    int(d['all_r_maps']), 
                    str(d['find_first_4leaves']), 
                    str(d['select_first_4leaves']),
                    int(d['common_r-steps']), 
                    int(d['common_r-steps_order']),
                ]

                if wp3_permu:
                    b = str(b).replace(',', '|')
                    res = res + [int(n_permu), int(max_permu), b]
                
                if small_spike:
                    res = res + [str(cycle)]
                
                res = res + [np.round(end - start, decimals=4)]
                df.loc[c, :] = res

                if c%l == 0 and c != 0:
                    with open(target, 'a') as f:
                        if wp3_permu:
                            seeds = f"# seed_szenario:{seed_szenario}|"\
                                    f" seed_pipeline:{seed_pipeline}|" \
                                    f" seed_permutations:{seed_permutations}\n"
                        else:
                            seeds = f"# seed_szenario:{seed_szenario}|"\
                                    f" seed_pipeline:{seed_pipeline}\n"

                        f.write(f"# Original Name: {file_name}\n")
                        f.write(seeds)
                        df.to_csv(f)
                    logger.info(f"Step: {c} / {r*N} - write DataFrame (partially)")
                    # df.to_csv(target)


                # if c == 8:
                #     name = f"extract_{file_name}_k={k}_i={i+1}"
                #     scenario.write_history(os.path.join(result_dir, name))

                c += 1
        
        # write final file
        if os.path.exists(target):
            os.remove(target)

        with open(target, 'a') as f:
            if wp3_permu:
                seeds = f"# seed_szenario:{seed_szenario}|"\
                    f" seed_pipeline:{seed_pipeline}|" \
                    f" seed_permutations:{seed_permutations}\n"
            else:
                seeds = f"# seed_szenario:{seed_szenario}|"\
                    f" seed_pipeline:{seed_pipeline}\n"

            f.write(f"# Original Name: {file_name}\n")
            f.write(seeds)
            df.to_csv(f)

        end_global = time.time()
        total_time_hour = (end_global - start_global) / 3600
        logger.info(f"Process completed. (it takes {np.round(total_time_hour, decimals=4)} h)")

    # ===========================================    
    def wp1(self, N:int=10, min_leafes:int=6, max_leaves:int=8, first_leaves:list=[0,1,2,3], 
        seed_szenario:int=0, add_name:str="") -> None:
        """ Execute workpackage 1 (See 'praktikumsanleitung_v.1.5').

        Parameter
        ---------
        N : int
            Number of generated random simulation szenarios with 'erdbeermet.simulation'. The
            options (circular, clocklike, first_candidate_only) must be set with the class
            'Pipeline' otherwise the default values are taken (circular=False, clocklike=False,
            first_candidate_only=True). The algorithmen will be executed for N simulations and
            all pissible number of leafes between min_leafes and max_leafes (both included).
            The default number of random simulations per number of leafes is 10.
        min_leafes : int, optional
            Minimal number of leaves for N simulations. The default is 6.
        max_leafes : int, optional
            Maximum number of leaves for N simulations. The default is 8.
        first_leaves : list, optional
            List of the first 4 leafes in the simulation. The default is [0,1,2,3] because
            the simulated scenarios from 'erdbeermet.simulation' will always start with 0 to 4.
        seed_szenario : int, optional
            Seed to make all szenarios for the same number of N simulations reproducible. By default
            it is set to 0 (if None, use a random seed between 0 and max int32).
        add_name : str, optional
            Additional name for the simulation process.

        Generates
        ---------
        wp1_results : folder
            Contains a csv file with all results for the simulations process (ex. if the recognition
            algorithmen run into an dead end).
        
        """
        if seed_szenario is None:
            seed_szenario = np.random.randint(np.iinfo(np.int32).max)

        self._workpackages(simulation="wp1",
            N=N, min_leaves=min_leafes, max_leaves=max_leaves, first_leaves=first_leaves,
            B=[], choose_random=False, small_spike=False, wp3_permu=False, subset=0,
            seed_szenario=seed_szenario, seed_pipeline=0,
            seed_permutations=0, add_name=add_name
        )
    

    # ===========================================
    def wp2(self, N:int=10, min_leafes:int=6, max_leaves:int=8, first_leaves:list=[0,1,2,3], 
        choose_random:bool=True, seed_szenario:int=0, seed_pipeline:int=0, add_name:str="") -> None:
        """ Execute workpackage 2 (See 'praktikumsanleitung_v.1.5').

        Parameter
        ---------
        N : int
            Number of generated random simulation szenarios with 'erdbeermet.simulation'. The
            options (circular, clocklike, first_candidate_only) must be set with the class
            'Pipeline' otherwise the default values are taken (circular=False, clocklike=False,
            first_candidate_only=True). The algorithmen will be executed for N simulations and
            all pissible number of leafes between min_leafes and max_leafes (both included).
            The default number of random simulations per number of leafes is 10.
        min_leafes : int, optional
            Minimal number of leaves for N simulations. The default is 6.
        max_leafes : int, optional
            Maximum number of leaves for N simulations. The default is 8.
        first_leaves : list, optional
            List of the first 4 leafes in the simulation. The default is [0,1,2,3] because
            the simulated scenarios from 'erdbeermet.simulation' will always start with leaves 
            identifier from 0 to 3. List will use to identify if the first four leaves are 
            correct recognized.
        choose_random : bool, optional
            If set to True, in the pipeline algorithm out of the list of candidates for 
            the last step, only one random candidate (with 4 leafes) is output as a positive
            R-Metric. All other will be set to 'random loser'. The default is True (if set
            to False it will produce the same output like 'wp1'.)
        seed_szenario : int, optional
            Seed to make all szenarios for the same number of N simulations reproducible. By default
            it is set to 0 (if None, use a random seed between 0 and max int32). 
        seed_pipeline : int, optional
            Seed to make all randomly choosen positive candidate for the same number of N simulations
            (with same seed_szenario) reproducible. Only useful with 'choose_random=True'. By 
            default it is set to 0 (if None, use a random seed between 0 and max int32).
        add_name : str, optional
            Additional name for the simulation process.

        Generates
        ---------
        wp2_results : folder
            Contains a csv file with all results for the simulations process (ex. if the recognition
            algorithmen run into an dead end).
        """
        if seed_szenario is None:
            seed_szenario = np.random.randint(np.iinfo(np.int32).max)

        if seed_pipeline is None:
            seed_pipeline = np.random.randint(np.iinfo(np.int32).max)

        self._workpackages(simulation="wp2",
            N=N, min_leaves=min_leafes, max_leaves=max_leaves, first_leaves=first_leaves,
            B=[], choose_random=choose_random, small_spike=False, wp3_permu=False, subset=4, 
            seed_szenario=seed_szenario, seed_pipeline=seed_pipeline, 
            seed_permutations=0, add_name=add_name
        )
    

    # ===========================================
    def wp3(self, N:int=10, min_leafes:int=6, max_leaves:int=8, first_leaves:list=[0,1,2,3], B:list=[0,1,2,3],
        leaves_unknown=False, subset:int=4, choose_random:bool=False, seed_szenario:int=0, seed_pipeline:int=0, 
        seed_permutations:int=0, add_name:str="") -> None:
        """ Execute workpackage 3 (See 'praktikumsanleitung_v.1.5').

        Parameter
        ---------
        N : int
            Number of generated random simulation szenarios with 'erdbeermet.simulation'. The
            options (circular, clocklike, first_candidate_only) must be set with the class
            'Pipeline' otherwise the default values are taken (circular=False, clocklike=False,
            first_candidate_only=True). The algorithmen will be executed for N simulations and
            all pissible number of leafes between min_leafes and max_leafes (both included).
            The default number of random simulations per number of leafes is 10.
        min_leafes : int, optional
            Minimal number of leaves for N simulations. The default is 6.
        max_leafes : int, optional
            Maximum number of leaves for N simulations. The default is 8.
        first_leaves : list, optional
            List of the first 4 leafes in the simulation. The default is [0,1,2,3] because
            the simulated scenarios from 'erdbeermet.simulation' will always start with leaves 
            identifier from 0 to 3.
        B : list, optional
            Leaf-identifier for the recognition algorithm. During the recognition no leaf
            from this list will be chosen to reduce as z. The default is [0,1,2,3] for this
            workpackage.
        leaves_unknown : bool, optional
            If set to True algorithm will iterate over permutations (in random order) in range of 
            k-leaves (k is equal to the current number of items) as leaf identifier B of size 
            'subset' (3 or 4) until an R-Map is correctly identified (reconition tree has no dead 
            end). In this case the result .csv file will have the following additional columns:
            'permu_count' (number of iteratios needed) 'max_permu' (maximum number of possible
            permutations) and 'permutation' (permutation as leaf_identifier B which produces the
            first R-Maps). The default is False.
        subset : int, optional
            Number of core leaves from which random permutations are created and used as 
            leav-identifier. Only useful if 'leaves_not_known=True'. Default is 4.
        choose_random : bool, optional
            If set to True, in the pipeline algorithm out of the list of candidates for 
            the last step, only one random candidate (with 4 leafes) is output as a positive
            R-Metric. All other will be set to 'random loser'. The default is False for this
            workpackage.
        seed_szenario : int, optional
            Seed to make all szenarios for the same number of N simulations reproducible. By default
            it is set to 0 (if None, use a random seed between 0 and max int32). 
        seed_pipeline : int, optional
            Seed to make all randomly choosen positive candidate for the same number of N simulations
            (with same seed_szenario) reproducible. Only useful with 'choose_random=True'. By 
            default it is set to 0 (if None, use a random seed between 0 and max int32).
        seed_permutations: int, optional
            Seed to make the order of permutations reproducible. Only useful if 'leaves_unknown=True'. 
            By default it is set to 0 (if None, use a random seed between 0 and max int32). 
        add_name : str, optional
            Additional name for the simulation process.

        Generates
        ---------
        wp3_results : folder
            Contains a csv file with all results for the simulations process (ex. if the recognition
            algorithmen run into an dead end).
        """



        if seed_szenario is None:
            seed_szenario = np.random.randint(np.iinfo(np.int32).max)

        if seed_pipeline is None:
            seed_pipeline = np.random.randint(np.iinfo(np.int32).max)
        
        if seed_permutations is None:
            seed_permutations = np.random.randint(np.iinfo(np.int32).max)

        if not leaves_unknown:
            simulation = "wp3"
        else:
            simulation = "wp3_permutations"

        self._workpackages(simulation=simulation,
            N=N, min_leaves=min_leafes, max_leaves=max_leaves, first_leaves=first_leaves,
            B=B, choose_random=choose_random, small_spike=False, wp3_permu=leaves_unknown, subset=subset, 
            seed_szenario=seed_szenario, seed_pipeline=seed_pipeline, seed_permutations=seed_permutations, 
            add_name=add_name
        )


    # ===========================================
    def wp4(self, N:int=10, min_leafes:int=6, max_leaves:int=8, first_leaves:list=[0,1,2,3], small_spike:bool=True, 
        choose_random:bool=False, seed_szenario:int=0, seed_pipeline:int=0, add_name:str="") -> None:
        """ Execute workpackage 4 (See 'praktikumsanleitung_v.1.5').

        Parameter
        ---------
        N : int
            Number of generated random simulation szenarios with 'erdbeermet.simulation'. The
            options (circular, clocklike, first_candidate_only) must be set with the class
            'Pipeline' otherwise the default values are taken (circular=False, clocklike=False,
            first_candidate_only=True). The algorithmen will be executed for N simulations and
            all pissible number of leafes between min_leafes and max_leafes (both included).
            The default number of random simulations per number of leafes is 10.
        min_leafes : int, optional
            Minimal number of leaves for N simulations. The default is 6.
        max_leafes : int, optional
            Maximum number of leaves for N simulations. The default is 8.
        first_leaves : list, optional
            List of the first 4 leafes in the simulation. The default is [0,1,2,3] because
            the simulated scenarios from 'erdbeermet.simulation' will always start with leaves 
            identifier from 0 to 3.
        small_spike : bool, optional
            If set to True, the recognition algorithmen will try to calculate the candidates for 
            R-Steps with respect to the smallest spike length (for remaining items > 5). If multiple 
            such candidates exist, the recognition will use an candidate at random. Additionally the
            algorithmen will add "_cycle-[True|False]" to dead ends (if they exists).
        choose_random : bool, optional
            If set to True, in the pipeline algorithm out of the list of candidates for 
            the last step, only one random candidate (with 4 leafes) is output as a positive
            R-Metric. All other will be set to 'random loser'. The default is False for this
            workpackage.
        seed_szenario : int, optional
            Seed to make all szenarios for the same number of N simulations reproducible. By default
            it is set to 0 (if None, use a random seed between 0 and max int32). 
        seed_pipeline : int, optional
            Seed to make all randomly choosen positive candidate for the same number of N simulations
            (with same seed_szenario) reproducible. Only useful with 'choose_random=True'. By 
            default it is set to 0 (if None, use a random seed between 0 and max int32).
        add_name : str, optional
            Additional name for the simulation process.

        Generates
        ---------
        wp4_results : folder
            Contains a csv file with all results for the simulations process (ex. if the recognition
            algorithmen run into an dead end).
        """
        if seed_szenario is None:
            seed_szenario = np.random.randint(np.iinfo(np.int32).max)

        if seed_pipeline is None:
            seed_pipeline = np.random.randint(np.iinfo(np.int32).max)

        self._workpackages(simulation="wp4",
            N=N, min_leaves=min_leafes, max_leaves=max_leaves, first_leaves=first_leaves,
            B=[], choose_random=choose_random, small_spike=small_spike, wp3_permu=False, subset=4, 
            seed_szenario=seed_szenario, seed_pipeline=seed_pipeline, 
            seed_permutations=0, add_name=add_name
        )
    