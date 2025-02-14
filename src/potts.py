import random
import numpy as np
from src.coupling import read_params
from src.util import seq2onehot, onehot2seq, int2onehot, onehot2int

class PottsModel:
    def __init__(self, J, h, spins, is_dna, is_gapped):
        """
        J: (num_nodes, num_nodes) numpy array of coupling strengths
        h: (num_nodes, num_spin_states) numpy array of external fields
        spins: (num_nodes, num_spin_states) numpy array of spin values
        """
        self.J = J
        self.h = h
        self.spins = spins
        self.is_dna = is_dna
        self.is_gapped = is_gapped
        self.num_nodes = J.shape[0]
        self.num_states = h.shape[1]

    def compute_energy(self):
        """
        compute energy of the current state (spins)
        """

        return - np.einsum('ij,ij->', self.h, self.spins) - np.einsum('ik,ijkl,jl->', self.spins, self.J, self.spins)

    def eval_energy(self, spins):
        """
        compute energy of the given state (spins)
        """
        return - np.einsum('ij,ij->', self.h, spins) - np.einsum('ik,ijkl,jl->', spins, self.J, spins)

    def compute_delta_energy(self, mutations):
        """
        mutations: list of tuple (int_wt, pos, int_mut)
        compute delta energy by mutations
        """
        delta_energy = 0
        for int_wt, pos, int_mut in mutations:
            delta_energy += self.h[pos, int_wt] - self.h[pos, int_mut]
            for i in range(self.num_nodes):
                if i != pos:
                    nuc_i = onehot2int(self.spins[i])
                    delta_energy += self.J[pos, i, nuc_i, int_wt] - self.J[pos, i, nuc_i, int_mut]
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
            spins = seq2onehot(params["target_seq"], max_len=len(params["target_seq"]), is_dna=is_dna, is_gapped = is_gapped),
            is_dna = is_dna, 
            is_gapped = is_gapped
            )


    
if __name__ == "__main__":
    model = PottsModel.build_from_file("data/params.txt")
    model.sim_anneal(1, 0.1, 1000)
    print(onehot2seq(model.spins))