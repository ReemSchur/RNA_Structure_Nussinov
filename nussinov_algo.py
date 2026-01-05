import numpy as np

class NussinovRNA:
    def __init__(self, sequence):
        """
        Initialize the NussinovRNA solver.
        """
        self.sequence = sequence.strip().upper().replace('T', 'U')
        self.n = len(self.sequence)
        self.dp_matrix = np.zeros((self.n, self.n), dtype=int)
        self.structure = ["."] * self.n

    def is_pair(self, base1, base2):
        pairs = {('A', 'U'), ('U', 'A'), 
                 ('G', 'C'), ('C', 'G'), 
                 ('G', 'U'), ('U', 'G')}
        return (base1, base2) in pairs

    def fill_matrix(self):
        for l in range(1, self.n):
            for i in range(self.n - l):
                j = i + l
                opt1 = self.dp_matrix[i + 1][j]
                opt2 = self.dp_matrix[i][j - 1]
                opt3 = 0
                if self.is_pair(self.sequence[i], self.sequence[j]):
                    opt3 = self.dp_matrix[i + 1][j - 1] + 1
                opt4 = 0
                for k in range(i, j):
                    current_split = self.dp_matrix[i][k] + self.dp_matrix[k + 1][j]
                    if current_split > opt4:
                        opt4 = current_split
                self.dp_matrix[i][j] = max(opt1, opt2, opt3, opt4)

    def traceback(self, i, j):
        if i >= j:
            return
        if self.dp_matrix[i][j] == self.dp_matrix[i + 1][j]:
            self.traceback(i + 1, j)
        elif self.dp_matrix[i][j] == self.dp_matrix[i][j - 1]:
            self.traceback(i, j - 1)
        elif (self.is_pair(self.sequence[i], self.sequence[j]) and 
              self.dp_matrix[i][j] == self.dp_matrix[i + 1][j - 1] + 1):
            self.structure[i] = "("
            self.structure[j] = ")"
            self.traceback(i + 1, j - 1)
        else:
            for k in range(i, j):
                if self.dp_matrix[i][j] == self.dp_matrix[i][k] + self.dp_matrix[k + 1][j]:
                    self.traceback(i, k)
                    self.traceback(k + 1, j)
                    break

    def predict(self):
        self.fill_matrix()
        self.structure = ["."] * self.n
        self.traceback(0, self.n - 1)
        return "".join(self.structure)

    def get_pair_count(self):
        """
        Returns the number of base pairs found in the predicted structure.
        """
        # Count the number of opening brackets '(' which represents the number of pairs
        return self.structure.count("(")
