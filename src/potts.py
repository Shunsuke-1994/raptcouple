import random
import numpy as np
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from plmc import read_params
from util import seq2onehot, onehot2seq, int2onehot, onehot2int

class PottsModel:
    def __init__(self, J, h, spins):
        """
        J: (num_nodes, num_nodes, num_state, num_state), assume symmetric and diagonal elements are zero
        h: (num_nodes, num_spin_states) numpy array of external fields
        spins: (num_nodes, num_spin_states) numpy array of spin values
        """
        self.J = J
        self.h = h
        self.spins = spins
        # self.is_dna = is_dna
        # self.is_gapped = is_gapped
        self.num_nodes = J.shape[0]
        self.num_states = h.shape[1]

    def compute_energy(self):
        """
        compute energy of the current state (spins)
        """

        return - np.einsum('ij,ij->', self.h, self.spins) - 0.5 * np.einsum('ik,ijkl,jl->', self.spins, self.J, self.spins)

    def eval_energy(self, spins):
        """
        compute energy of the given state (spins)
        """
        return - np.einsum('ij,ij->', self.h, spins) - 0.5 * np.einsum('ik,ijkl,jl->', spins, self.J, spins)

    def compute_delta_energy(self, mutations, print_debug = False):
        """
        mutations: list of tuple (int_wt, pos, int_mut)
        compute delta energy by mutations
        """
        delta_energy = 0
        for int_wt, pos, int_mut in mutations:
            delta_position = self.h[pos, int_wt] - self.h[pos, int_mut]
            delta_energy += delta_position
            if print_debug:
                print(f"delta_position: {delta_position}, pos = {pos}, int_wt = {int_wt}, int_mut = {int_mut}")
            for i in range(self.num_nodes):
                if i != pos:
                    nuc_i = onehot2int(self.spins[i])

                    # from symmetry of J, we assume i < pos and double the energy(0.5x2)
                    delta_coupling = self.J[pos, i, int_wt, nuc_i] - self.J[pos, i, int_mut, nuc_i]
                    delta_energy += delta_coupling
                    if print_debug:
                        print(f"delta_coupling: {delta_coupling}, pos = {pos}, i = {i}, nuc_i = {nuc_i}, int_wt = {int_wt}, int_mut = {int_mut}")
        return delta_energy

    def sim_anneal(self, T_init, T_end, max_steps, random_init_state = False, random_seed = 42):
        """
        Simutated annealing by Metroplis-Hastings algorithm.
        T_init: initial temperature
        T_end: final temperature
        random_init_state: whether to randomly initialize the state
        max_steps: maximum number of steps
        """
        random.seed(random_seed)
        np.random.seed(random_seed)
        T = T_init
        if random_init_state:
            self.spins = np.zeros(self.spins.shape)
            spin_init = np.random.randint(self.num_states, size=self.num_nodes)
            for i in range(self.num_nodes):
                self.spins[i] = int2onehot(spin_init[i], self.num_states)
        current_energy = self.compute_energy()

        for step in range(max_steps):
            # randomly choose a node to flip
            node = np.random.randint(self.num_nodes)
            current_spin_int = onehot2int(self.spins[node])
            # sample other than current_spin states
            new_spin_int = np.random.choice([i for i in range(self.num_states) if i != current_spin_int])

            # compute delta energy
            delta_energy = self.compute_delta_energy([(current_spin_int, node, new_spin_int)])

            # lower energy -> accept 
            if delta_energy < 0: # accept
                self.spins[node] = int2onehot(new_spin_int, self.num_states)
                current_energy += delta_energy
            else: # accept with probability
                if np.random.rand() < np.exp(-delta_energy / T): # r ~ U[0,1]
                    self.spins[node] = int2onehot(new_spin_int, self.num_states)
                    current_energy += delta_energy
            # update temperature
            T = T_init * (T_end / T_init) ** (step / max_steps)

        return 
    
    def gibbs_sampling(self, n_samples, random_seed = 42):
        """
        Gibbs sampling
        n_samples: number of samples
        """
        random.seed(random_seed)
        np.random.seed(random_seed)
        
        samples = np.zeros((n_samples, self.num_nodes, self.num_states))
        for n in range(n_samples):
            for i in range(self.num_nodes):
                # compute probabilities of each state
                probs = np.exp([
                    -self.compute_delta_energy([(onehot2int(self.spins[i]), i, k)])
                    for k in range(self.num_states)
                    ])
                probs /= np.sum(probs)
                # sample a state
                samples[n, i] = np.random.multinomial(1, probs)
        return samples

    @classmethod
    def build_from_file(cls, param_file):
        params = read_params(param_file)
        is_dna = True if "T" in params["target_seq"] else False
        is_gapped = True if "." in params["alphabet"] else False

        return cls(
            J = params["Jij"],
            h = params["hi"],
            spins = seq2onehot(params["target_seq"], max_len=len(params["target_seq"]), is_dna=is_dna, is_gapped = is_gapped)
            )


    
if __name__ == "__main__":
    # Original test code
    n_nodes = 3    # ノード数
    n_states = 2   # 各ノードの状態数

    # J: (shape: (n_nodes, n_nodes, n_states, n_states))
    # random J, synmetric, diagonal elements are zero
    J = np.random.rand(n_nodes, n_nodes, n_states, n_states)
    for i in range(n_nodes):
        for j in range(i+1, n_nodes):
            # J[j, i] = (J[i, j])の転置で対称性を確保
            J[j, i] = J[i, j].T
    # 対角成分はゼロに設定
    for i in range(n_nodes):
        J[i, i] = np.zeros((n_states, n_states))

    # 外部場 h の作成 (shape: (n_nodes, n_states))
    h = np.random.rand(n_nodes, n_states)
    print("J:", J)
    print("h:", h)
    # 初期スピン状態の定義 (例: ノードごとに状態 [1, 0] または [0, 1])
    spins = np.array([
        [1, 0],
        [0, 1],
        [1, 0]
        ])

    model = PottsModel(J, h, spins)
    energy_direct = model.compute_energy()
    energy_eval = model.eval_energy(model.spins)
    energy_delta_null = model.compute_delta_energy([(0, 0, 0)])
    print("Test eval_energy:")
    print("Energy (compute_energy):", energy_direct)
    print("Energy (eval_energy):", energy_eval)
    print("Null delta:", energy_delta_null)  # check if the two methods give the same result
    print("Difference:", energy_eval - energy_direct)  # check if the two methods give the same result

    # --- compute_delta_energy のテスト ---
    # ここでは、あるノードのスピンを1つ変化させた場合のエネルギー変化が、
    # eval_energy を用いて直接計算したエネルギー差と一致するかを確認する。
    print("\nTest compute_delta_energy:")
    pos = 1  # 例としてノード1を選択
    current_spin = onehot2int(model.spins[pos])
    new_spin = 1 - current_spin

    # deltaE by compute_delta_energy 
    delta_energy_computed = model.compute_delta_energy([(current_spin, pos, new_spin)], print_debug = True)

    new_spins = model.spins.copy()
    new_spins[pos] = int2onehot(new_spin, n_states)
    print(model.spins[pos], "->", new_spins[pos], "at node", pos)
    print(model.spins, "->", new_spins)
    energy_new = model.eval_energy(new_spins)
    energy_old = energy_direct

    # deltaE by eval_energy
    delta_energy_actual = energy_new - energy_old
    print("Energy before:", energy_old)
    print("Energy after :", energy_new)

    print("deltaE by compute_delta_energy:\t", delta_energy_computed)
    print("deltaE by eval_energy:\t\t", delta_energy_actual)
    print("Difference of compute_delta_energy and eval_energy:", delta_energy_computed - delta_energy_actual)


