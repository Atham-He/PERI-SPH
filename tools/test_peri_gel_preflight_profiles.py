import unittest
from peri_gel_preflight import validate_limits


class GelPreflightProfiles(unittest.TestCase):
    def test_legacy_profile_keeps_original_1000_bound(self):
        validate_limits(dict(max_particles=1000, default_threads=8))

    def test_existing_spatial_profile_requires_its_explicit_identity(self):
        validate_limits(dict(execution_profile='solid-particle-spatial-v1', max_particles=8192,
                             max_steps=98304, max_threads=8, default_threads=8))

    def test_unlabeled_or_unknown_enlargement_rejects(self):
        for extra in [{}, {'execution_profile': 'unknown'}]:
            with self.assertRaises(ValueError):
                validate_limits(dict(max_particles=8192, default_threads=8, **extra))

    def test_modified_declared_limits_and_bad_thread_counts_reject(self):
        base = dict(execution_profile='solid-particle-spatial-v1', max_particles=8192,
                    max_steps=98304, max_threads=8, default_threads=8)
        for key, value in [('max_particles', 8193), ('max_steps', 98305), ('max_threads', 9),
                           ('default_threads', 9), ('default_threads', 0)]:
            with self.assertRaises(ValueError):
                validate_limits(dict(base, **{key: value}))


if __name__ == '__main__':
    unittest.main()
