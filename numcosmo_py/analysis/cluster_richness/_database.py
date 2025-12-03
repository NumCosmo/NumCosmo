#
# _database.py
#
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""BestfitDatabase for ACID-compliant storage of analysis results.

This module provides SQLite-based storage with crash safety guarantees
for best-fit parameters from cluster mass-richness analysis.

The database now uses NumCosmo's serialization for storing complete model
instances, allowing any NcClusterMassRichness subclass to be stored.
"""

import sqlite3

from numcosmo_py import Nc

from ._parameters import (
    model_params_to_dict,
    model_to_yaml,
    model_from_yaml,
    CutAnalysisResult,
)


class BestfitDatabase:
    """SQLite database manager for best-fit models with ACID guarantees.

    Provides efficient storage and retrieval of best-fit models indexed by
    mock seed and cut, with automatic recovery on script restart and crash safety
    via transaction wrapping.

    The database stores serialized model objects, allowing any
    NcClusterMassRichness subclass to be stored and retrieved.
    """

    def __init__(self, db_path: str = "bestfits.db"):
        """Initialize database connection and create schema if needed.

        :param db_path: Path to SQLite database file (default: "bestfits.db")
        """
        self.db_path = db_path
        self._init_db()

    def _init_db(self):
        """Initialize database schema with ACID-compliant transaction handling."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                "PRAGMA journal_mode = WAL"
            )  # Write-Ahead Logging for crash safety
            conn.execute("PRAGMA synchronous = NORMAL")  # Balance speed and safety
            cursor = conn.cursor()

            # Create table with composite unique key (mock_seed, cut)
            # Stores serialized model as YAML string for flexibility
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS best_fits (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    mock_seed INTEGER NOT NULL,
                    cut REAL NOT NULL,
                    n_clusters INTEGER,
                    model_type TEXT NOT NULL,
                    model_yaml TEXT NOT NULL,
                    m2lnL REAL NOT NULL,
                    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
                    UNIQUE(mock_seed, cut)
                )
                """
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_mock_seed ON best_fits(mock_seed)"
            )
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_cut ON best_fits(cut)")
            conn.commit()

    def insert_bestfit(self, mock_seed: int, result: CutAnalysisResult) -> bool:
        """Insert a best-fit result, with transaction wrapping for crash safety.

        :param mock_seed: The mock seed identifier
        :param result: CutAnalysisResult containing best-fit model
        :return: True if inserted, False if already exists
        """
        try:
            yaml_str = model_to_yaml(result.bestfit)
            model_type = type(result.bestfit).__name__

            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute(
                    """
                    INSERT OR IGNORE INTO best_fits
                    (mock_seed, cut, n_clusters, model_type, model_yaml, m2lnL)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        mock_seed,
                        result.cut,
                        result.n_clusters,
                        model_type,
                        yaml_str,
                        result.m2lnL,
                    ),
                )
                conn.commit()
                return cursor.rowcount > 0
        except sqlite3.IntegrityError:
            return False

    def batch_insert(self, seeds_results: list[tuple[int, CutAnalysisResult]]):
        """Batch insert multiple results in a single transaction for efficiency.

        :param seeds_results: List of (mock_seed, CutAnalysisResult) tuples
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            for mock_seed, result in seeds_results:
                yaml_str = model_to_yaml(result.bestfit)
                model_type = type(result.bestfit).__name__

                cursor.execute(
                    """
                    INSERT OR IGNORE INTO best_fits
                    (mock_seed, cut, n_clusters, model_type, model_yaml, m2lnL)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        mock_seed,
                        result.cut,
                        result.n_clusters,
                        model_type,
                        yaml_str,
                        result.m2lnL,
                    ),
                )
            conn.commit()

    def get_computed_seeds(self, cut: float) -> set[int]:
        """Get set of already-computed seeds for a specific cut.

        :param cut: The richness cut value
        :return: Set of mock_seed values already in database
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT mock_seed FROM best_fits WHERE cut = ?", (cut,))
            return {row[0] for row in cursor.fetchall()}

    def get_all_computed_seeds(self) -> set[int]:
        """Get set of all seeds with at least one cut computed.

        :return: Set of all mock_seed values in database
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT DISTINCT mock_seed FROM best_fits")
            return {row[0] for row in cursor.fetchall()}

    def get_bestfit(self, mock_seed: int, cut: float) -> Nc.ClusterMassRichness | None:
        """Retrieve best-fit model for a specific seed and cut.

        :param mock_seed: The mock seed identifier
        :param cut: The richness cut value
        :return: Deserialized model or None if not found
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            cursor.execute(
                "SELECT model_yaml FROM best_fits WHERE mock_seed = ? AND cut = ?",
                (mock_seed, cut),
            )
            row = cursor.fetchone()
            if row:
                return model_from_yaml(row["model_yaml"])
            return None

    def get_bestfit_params_dict(
        self, mock_seed: int, cut: float
    ) -> dict[str, float] | None:
        """Retrieve best-fit parameters as a dictionary.

        :param mock_seed: The mock seed identifier
        :param cut: The richness cut value
        :return: Parameter dictionary or None if not found
        """
        model = self.get_bestfit(mock_seed, cut)
        if model is not None:
            return model_params_to_dict(model)
        return None

    def delete_bestfit(self, mock_seed: int, cut: float) -> bool:
        """Delete a best-fit entry (for recomputation with --recompute flag).

        :param mock_seed: The mock seed identifier
        :param cut: The richness cut value
        :return: True if deleted, False if not found
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(
                "DELETE FROM best_fits WHERE mock_seed = ? AND cut = ?",
                (mock_seed, cut),
            )
            conn.commit()
            return cursor.rowcount > 0

    def count_entries(self) -> int:
        """Get total number of best-fit entries in database.

        :return: Number of entries
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM best_fits")
            return cursor.fetchone()[0]
