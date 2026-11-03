import os
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import peri_run


class SourceSnapshotTests(unittest.TestCase):
    def test_build_profiles_and_authored_assets_are_frozen_but_runtime_trees_are_not(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            kept = ['CMakeLists.txt', 'src/shared/peridynamics/pd_fsi.h', 'src/shared/output/observer.h',
                    'peri/CMakeLists.txt', 'peri/experiments/paper.json', 'tools/peri_experiment.py',
                    'tests/3d_examples/extra_3d_PD_case/CMakeLists.txt',
                    'tests/3d_examples/extra_3d_PD_case/input/shape.stl',
                    'tests/3d_examples/extra_3d_PD_case/reload/particles.xml']
            omitted = ['src/shared/__pycache__/cache.pyc', 'peri/.venv/library.py',
                       'tests/3d_examples/extra_3d_PD_case/output/large_raw.csv',
                       'tests/3d_examples/extra_3d_PD_case/checkpoints/state.bin',
                       'tests/3d_examples/extra_3d_PD_case/build/program.cpp',
                       'tests/3d_examples/extra_3d_PD_case/generated.o']
            for name in kept + omitted:
                path = repo / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('source or generated data')
            actual = {p.relative_to(repo).as_posix() for p in peri_run.source_snapshot_files(repo)}
            self.assertEqual(actual, set(kept))

    def test_generated_runtime_directory_is_pruned_before_descent(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            case = repo / 'tests/3d_examples/extra_3d_PD_case'
            (case / 'output/deep').mkdir(parents=True)
            (case / 'output/deep/large.csv').write_text('raw')
            visited = []
            original_walk = os.walk
            def tracked_walk(path, **kwargs):
                for item in original_walk(path, **kwargs):
                    visited.append(Path(item[0]))
                    yield item
            with patch('peri_run.os.walk', tracked_walk):
                peri_run.source_snapshot_files(repo)
            self.assertNotIn(case / 'output', visited)
            self.assertNotIn(case / 'output/deep', visited)

    def test_unreadable_source_directory_does_not_become_a_partial_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            (repo / 'src').mkdir()
            def unreadable(path, **kwargs):
                kwargs['onerror'](PermissionError('cannot read source'))
            with patch('peri_run.os.walk', unreadable):
                with self.assertRaisesRegex(PermissionError, 'cannot read source'):
                    peri_run.source_snapshot_files(repo)

    def test_stale_source_fails_before_creating_or_freezing_a_run(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            (repo / 'tools').mkdir()
            source = repo / 'src/shared/example.cpp'
            source.parent.mkdir(parents=True)
            source.write_text('new source')
            executable = repo / 'solver'
            executable.write_text('old executable; must never be invoked')
            config = repo / 'case.cfg'
            config.write_text('with_liquid=true\n')
            os.utime(executable, (10, 10))
            os.utime(source, (20, 20))
            run = repo / 'new_run'
            args = SimpleNamespace(executable=executable, config=config, run_dir=run)
            with patch.object(peri_run, '__file__', str(repo / 'tools/peri_run.py')):
                with self.assertRaisesRegex(ValueError, 'rebuild before running'):
                    peri_run.launch(args)
            self.assertFalse(run.exists())


if __name__ == '__main__':
    unittest.main()
