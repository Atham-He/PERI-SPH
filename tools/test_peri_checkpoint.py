import json
import os
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from peri_run import committed_checkpoints,digest,resume,solver_alive
from peri_control import request,restart_study
import test_peri_audit as audit_fixtures


class CheckpointControlTests(unittest.TestCase):
    def test_partial_or_truncated_generations_are_not_committed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name,size in [("checkpoint-1",6),("checkpoint-2",2),(".tmp-checkpoint-3",6)]:
                p=root/name;p.mkdir();(p/"state.bin").write_bytes(b"abcdef"[:size])
                (p/"manifest.json").write_text(json.dumps({"format":1,"status":"complete","state_bytes":6}))
            self.assertEqual([Path(x["path"]).name for x in committed_checkpoints(root)],["checkpoint-1"])

    @unittest.skipUnless(Path("/proc/sys/kernel/random/boot_id").exists(),"Linux process identity test")
    def test_pid_start_and_boot_must_match_before_control(self):
        fields=(Path("/proc")/str(os.getpid())/"stat").read_text().rsplit(")",1)[1].split()
        metadata={"solver_pid":os.getpid(),"linux_process_start_ticks":fields[19],"linux_boot_id":Path("/proc/sys/kernel/random/boot_id").read_text().strip()}
        self.assertTrue(solver_alive(metadata))
        self.assertFalse(solver_alive({**metadata,"linux_process_start_ticks":"0"}))
        self.assertFalse(solver_alive({**metadata,"linux_boot_id":"different-boot"}))

    def test_resume_rejects_changed_binary_without_executing(self):
        with tempfile.TemporaryDirectory() as directory:
            run=Path(directory);exe=run/"solver";exe.write_bytes(b"original");config=run/"case.cfg";config.write_text("parameter=1\n")
            metadata={"frozen_executable":str(exe),"executable_sha256":digest(exe),"case_sha256":digest(config),"checkpoint_format":1,"state_recording":False}
            (run/"run_manifest.json").write_text(json.dumps(metadata));exe.write_bytes(b"modified")
            with patch("peri_run.execute") as execute:
                with self.assertRaisesRegex(ValueError,"sealed executable"):resume(SimpleNamespace(run_dir=run,checkpoint=None))
                execute.assert_not_called()

    def test_legacy_pause_never_sends_an_unsupported_signal(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/"runs/baseline").mkdir(parents=True)
            (root/"study.json").write_text(json.dumps({"active_case":"baseline"}))
            (root/"runs/baseline/run_manifest.json").write_text(json.dumps({"checkpoint_format":0}))
            with patch("peri_control.os.kill") as kill:
                with self.assertRaisesRegex(ValueError,"no full PD checkpoint"):request(root,True,0)
                kill.assert_not_called()

    def test_disk_pause_is_not_completed_or_a_dead_solver_warning(self):
        result=audit_fixtures.AcousticBudgetTests().audit(state="checkpointed\ntime_s=1")
        self.assertTrue(result["simulation_paused_on_disk"])
        self.assertFalse(result["simulation_completed"])
        self.assertEqual(result["status"],"sampled_checks_passed")

    def test_missing_checkpoint_cannot_silently_restart_study_from_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/"runs/baseline").mkdir(parents=True)
            (root/"study.json").write_text(json.dumps({"status":"paused","active_case":"baseline"}))
            (root/"runs/baseline/run_manifest.json").write_text(json.dumps({"checkpoint_format":1}))
            with patch("peri_control.subprocess.Popen") as launch:
                with self.assertRaisesRegex(ValueError,"refusing a silent restart from zero"):restart_study(root,True)
                launch.assert_not_called()

if __name__=="__main__":unittest.main()
