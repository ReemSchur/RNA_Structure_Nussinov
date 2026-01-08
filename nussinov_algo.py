import numpy as np

class NussinovRNA:
    def __init__(self, sequence):
        """
        Initialize the NussinovRNA solver.
        """
        self.sequence = sequence.strip().upper().replace('T', 'U').replace(" ", "")
        self.n = len(self.sequence)
        
        # M(i,j) stores the max number of pairs in subsequence i..j
        self.dp_matrix = np.zeros((self.n, self.n), dtype=int)
        
        # K(i,j) stores the index 'k' that paired with 'j' to give the max score.
        # Initialize with -1 to represent "no pair for j".
        self.k_matrix = np.full((self.n, self.n), -1, dtype=int)
        
        self.structure = ["."] * self.n

    def is_pair(self, base1, base2):
        """
        Checks for Watson-Crick pairs (A-U, G-C).
        Adjacent pairs are allowed (no minimum distance check).
        """
        pairs = {('A', 'U'), ('U', 'A'), 
                 ('G', 'C'), ('C', 'G')}
        return (base1, base2) in pairs

    def fill_matrix(self):
        """
        Fills M(i,j) and K(i,j) according to Nussinov & Jacobson (1980).
        Logic follows Equation [1] in the paper.
        """
        # Iterate over subsequence lengths
        for l in range(1, self.n):
            for i in range(self.n - l):
                j = i + l
                
                # Option 1: j is unpaired (Skip j)
                # M(i,j) = M(i, j-1)
                max_score = self.dp_matrix[i][j-1]
                best_k = -1 # -1 indicates j is unpaired
                
                # Option 2: j is paired with some k in range [i, j-1]
                # Formula: M(i, k-1) + M(k+1, j-1) + 1
                for k in range(i, j):
                    if self.is_pair(self.sequence[k], self.sequence[j]):
                        # Calculate score for splitting at k
                        # Careful with indices: valid ranges only
                        left_part = self.dp_matrix[i][k-1] if k > i else 0
                        inner_part = self.dp_matrix[k+1][j-1] if k+1 <= j-1 else 0
                        
                        current_score = left_part + inner_part + 1
                        
                        if current_score > max_score:
                            max_score = current_score
                            best_k = k
                
                self.dp_matrix[i][j] = max_score
                self.k_matrix[i][j] = best_k

    def traceback(self, i, j):
        """
        Reconstructs structure using the K matrix.
        See paper section: "This information is used to find the best folded form..."
        """
        if i >= j:
            return

        k = self.k_matrix[i][j]

        if k == -1:
            # Case: j is unpaired.
            # "If B_j cannot pair... then M(i,j) = M(i, j-1)"
            # Simply move j to the left.
            self.traceback(i, j - 1)
        else:
            # Case: j is paired with k.
            # "The first base pair formed is B_n B_K(1,n)" (adapted to current i,j)
            self.structure[k] = "("
            self.structure[j] = ")"
            
            # "Divides the sequence into two subsections... B_1...B_k-1 and B_k+1...B_n-1"
            # Recurse on the left of k (if any)
            if k > i:
                self.traceback(i, k - 1)
            
            # Recurse inside the pair (between k and j)
            if k + 1 < j:
                self.traceback(k + 1, j - 1)

    def predict(self):
        """
        Main execution.
        """
        self.fill_matrix()
        self.structure = ["."] * self.n
        # Start traceback from the full sequence (0 to n-1)
        self.traceback(0, self.n - 1)
        return "".join(self.structure)

    def get_pair_count(self):
        return self.structure.count("(")