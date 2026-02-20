from qiskit.quantum_info import SparsePauliOp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import qctrlvisualizer as qv
from scipy.sparse.linalg import eigsh


def clean(arr, tol=1e-15):
    arr = arr.copy()
    arr[np.abs(arr) < tol] = 0.0
    return arr


N = 12 # qubit number

# --- Lattice structure as a graph ---
def draw_lattice():
    G = nx.Graph()
    edges = [(i, (i + 1) % N) for i in range(N)]
    G.add_edges_from(edges)
    pos = {i: np.array([np.cos(-2 * np.pi * i / N + np.pi/2), 
                        np.sin(-2 * np.pi * i / N + np.pi/2)]) for i in range(N)}

    plt.figure(figsize=(5, 5))
    nx.draw(G, pos, with_labels=True, node_color=qv.QCTRL_STYLE_COLORS[0], node_size=800)
    plt.title(f"TFIM Lattice N={N}")
    plt.show()

# --- Lowest two energy levels computation ---

#Building the hamiltonian
def get_hamiltonian(N, g, sparse=True):
    
    sparse_list = []
    for i in range(N):
        sparse_list.append(("ZZ", [i, (i + 1) % N], -1.0))
    for i in range(N):
        sparse_list.append(("X", [i], -g))

    hamiltonian = SparsePauliOp.from_sparse_list(sparse_list, num_qubits=N)
    return hamiltonian.to_matrix(sparse=sparse)

#Getting the lowest energy levels
def compute_lowest_energies(N, g, tol=1e-15):

    H_sparse = get_hamiltonian(N, g, sparse=True)
    
    #Eigenvalues for the energies and ground state eigenvector for the structure factor S
    vals, vecs = eigsh(H_sparse, k=2, which='SA', return_eigenvectors=True)
    order = np.argsort(vals)
    psi0 = vecs[:, order[0]]
    psi0[np.abs(psi0) < tol] = 0.0
    return vals[order], psi0


def print_hamiltonian_matrix(N, g):
    #Visualizing the hamiltonian matrix cuz why not
    H_mat = np.real(get_hamiltonian(N, g, sparse=False))
    
    np.set_printoptions(precision=3, suppress=True, linewidth=250, threshold=np.inf)
    print(f"\nHamiltonian Matrix (N={N}, g={g}):")
    print(H_mat)

## -- structure factor ---
def get_szz_factor(N, g):

    ''' Since we're calling compute_lowest_energies, 
        might as well store the energies so that we don't 
        call it again for plotting them'''
        
    vals, psi0 = compute_lowest_energies(N, g)
    terms = []
    for i in range(N):
        for j in range(N):
            if i == j:
                pass  # Zi*Zi=I
            else:
                terms.append(("ZZ", [i, j], 1.0))
    
    ZZ_total = SparsePauliOp.from_sparse_list(terms, num_qubits=N).to_matrix(sparse=True)
    
    # N addition contributions Zi*Zi
    S = np.real(np.vdot(psi0, ZZ_total @ psi0)) + N
    
    return S / N**2, vals[0], vals[1]


if __name__ == "__main__":
    
    # --- Parameters for the sweep ---
    g_values = np.linspace(0, 3, 50)
    N_list = [4, 6, 8, 10, 12, 14]
    
    fig_S, ax_S = plt.subplots(figsize=(8, 6))
    fig_E, ax_E = plt.subplots(figsize=(8, 6))
    
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(N_list)))

    for i, N in enumerate(N_list):
        S_values = np.zeros_like(g_values)
        E0_values = np.zeros_like(g_values)
        E1_values = np.zeros_like(g_values)
        
        for idx, g in enumerate(g_values):
            S, E0, E1 = get_szz_factor(N, g)
            S_values[idx] = S
            E0_values[idx] = E0 / N
            E1_values[idx] = E1 / N
            
        ax_S.plot(g_values, S_values, color=colors[i], linewidth=2.2, label=f"N = {N}")
        
        if N == 6:
            ax_E.plot(g_values, E0_values, color='b', linestyle='-', linewidth=2.2, label=f"E0/N (N={N})")
            ax_E.plot(g_values, E1_values, color='g', linestyle='-', linewidth=2.2, label=f"E1/N (N={N})")

    ax_S.axvline(x=1.0, color="#FF5722", linestyle="--", alpha=0.7, label=r"$g_c = 1$")
    ax_S.set_xlabel(r"Transverse field strength $g$", fontsize=13)
    ax_S.set_ylabel(r"Structure factor $S$", fontsize=13)
    ax_S.set_title("Magnetic structure factor vs g for diff N", fontsize=14)
    ax_S.legend(fontsize=11)
    ax_S.grid(True, alpha=0.3)

    ax_E.axvline(x=1.0, color="#FF5722", linestyle="--", alpha=0.7, label=r"$g_c = 1$")
    ax_E.set_xlabel(r"Transverse field strength $g$", fontsize=13)
    ax_E.set_ylabel(r"Energy Density $E/N$", fontsize=13)
    ax_E.set_title("Lowest energies density vs g (N=6)", fontsize=14)
    ax_E.legend(fontsize=9, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax_E.grid(True, alpha=0.3)

    plt.tight_layout()  
    
    fig_S.savefig("figures/structure_factor.png", dpi=300, bbox_inches='tight')
    fig_E.savefig("figures/energy_density_N=6.png", dpi=300, bbox_inches='tight')
    
   
