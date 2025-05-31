# entraine.py
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression # Tel que décrit dans automatique.pdf [cite: 281]
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import joblib
import os

def train_model():
    model_filename = "ml_model.pkl"
    data_filename = "branching_data.csv"

    try:
        if not os.path.exists(data_filename):
            print(f"Erreur: {data_filename} non trouvé.")
            joblib.dump((None, None), model_filename); print(f"Fichier {model_filename} placeholder créé."); return False

        df = pd.read_csv(data_filename)
        df.columns = df.columns.str.strip()
        print(f"Données chargées pour entraînement: {len(df)} échantillons")

        features_for_X = ['y_value', 'node_depth']
        target_col = 'score_strong_branching' # La cible est le score du SB

        required_cols = features_for_X + [target_col]
        for col in required_cols:
            if col not in df.columns:
                print(f"Erreur: Colonne '{col}' manquante dans {data_filename}.")
                joblib.dump((None, None), model_filename); print(f"Fichier {model_filename} placeholder créé."); return False
        
        for col in df.columns:
            if col in required_cols: # Seulement convertir les colonnes nécessaires
                 df[col] = pd.to_numeric(df[col], errors='coerce')

        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=required_cols)
        
        print(f"Données après nettoyage: {len(df)} échantillons")

        if len(df) < 10: # Seuil pour un entraînement minimal
            print("Pas assez de données valides pour l'entraînement (moins de 10).")
            joblib.dump((None, None), model_filename); print(f"Fichier {model_filename} placeholder créé."); return True

        X = df[features_for_X].values
        y = df[target_col].values

        if len(df) < 50: # Si peu de données, entraîner sur tout et normaliser sur tout
            print("Entraînement sur toutes les données (peu d'échantillons).")
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            model = LinearRegression()
            model.fit(X_scaled, y)
            score_val = model.score(X_scaled, y) 
            joblib.dump((model, scaler), model_filename)
            print(f"Modèle et scaler sauvegardés dans {model_filename}.")
            return True

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        model = LinearRegression()
        model.fit(X_train_scaled, y_train) # Apprentissage des coefficients beta [cite: 286]

        train_score_val = model.score(X_train_scaled, y_train)
        test_score_val = model.score(X_test_scaled, y_test)
       
        joblib.dump((model, scaler), model_filename)
        print(f"Modèle sauvegardé dans {model_filename}")
        return True

    except FileNotFoundError:
        print(f"Erreur: Fichier de données '{data_filename}' non trouvé."); joblib.dump((None, None), model_filename); print(f"Fichier {model_filename} placeholder créé."); return False
    except Exception as e:
        print(f"Erreur inattendue: {e}"); joblib.dump((None, None), model_filename); print(f"Fichier {model_filename} placeholder créé."); return False

if __name__ == "__main__":
    if train_model(): exit(0)
    else: exit(1)