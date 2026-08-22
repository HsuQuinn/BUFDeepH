from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class ProjectContractTest(unittest.TestCase):
    def test_band_unfolding_exposes_v310_cli_options(self):
        src = read("src/band_unfolding.jl")

        for option in [
            "--m-trans",
            "--num-band",
            "--max-iter",
            "--num-points",
            "--fermi-shift",
            "--point-start",
            "--point-end",
            "--part-index",
            "--num-parts",
            "--mid-z",
            "--output-file",
        ]:
            self.assertIn(option, src)

    def test_band_unfolding_uses_parameterized_latest_experiment_defaults(self):
        src = read("src/band_unfolding.jl")

        self.assertIn('BUFDeepH_VERSION = "3.1.0"', src)
        self.assertIn("M_TRANS_389", src)
        self.assertIn("M_TRANS_943", src)
        self.assertIn("partition_indices", src)
        self.assertIn("H_R, S_R, norbits, fermi_level, orbital_position = preprocess", src)
        self.assertNotIn("p_points = pp_points[45:49]", src)
        self.assertIn("function SolveHk(k, H_R, S_R, norbits, fermi_level", src)

    def test_slurm_array_script_runs_partitions_from_cli(self):
        script_path = ROOT / "scripts" / "slurm_band_unfolding_array.sh"
        self.assertTrue(script_path.exists())

        script = script_path.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --array=", script)
        self.assertIn("SLURM_ARRAY_TASK_ID", script)
        self.assertIn("--part-index", script)
        self.assertIn("--num-parts", script)
        self.assertIn("src/band_unfolding.jl", script)


if __name__ == "__main__":
    unittest.main()
