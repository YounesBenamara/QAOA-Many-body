from qiskit_ibm_runtime import QiskitRuntimeService

try:
    print("Tentative de connexion à IBM Quantum...")
    # Va chercher automatiquement le compte sauvegardé sur votre machine
    service = QiskitRuntimeService()
    
    # On récupère la liste des machines disponibles
    backends = service.backends()
    
    print(f"\n✅ Connexion réussie ! Vous avez accès à {len(backends)} backends Qiskit Runtime.")
    print("Voici quelques exemples de machines disponibles pour vous :")
    for b in backends[:3]:
        print(f" - {b.name}")
        
except Exception as e:
    print("\n❌ Erreur de connexion. Vos identifiants (token) sont peut-être invalides ou n'ont pas été correctement sauvegardés.")
    print(f"Détail de l'erreur : {e}")
