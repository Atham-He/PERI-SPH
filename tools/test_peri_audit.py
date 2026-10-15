from pathlib import Path
import tempfile
import unittest
from xml.parsers import expat

from peri_audit import rows, scan_vtp


VTP = '''<VTKFile><PolyData><FieldData><DataArray Name="TimeValue" type="Float64" NumberOfTuples="1" format="ascii">1.234567890123456e-11</DataArray></FieldData>
<Piece NumberOfPoints="1"><Points><DataArray Name="Position" type="Float64" NumberOfComponents="3" format="ascii">1e-9 2e-9 3e-9</DataArray></Points><PointData>
<DataArray Name="Mass" type="Float64" format="ascii">9.350000000000002e-13</DataArray>
<DataArray Name="Density" type="Float64" format="ascii">935</DataArray>
<DataArray Name="VolumetricMeasure" type="Float64" format="ascii">1.0000000000000003e-15</DataArray>
<DataArray Name="Velocity" type="Float64" NumberOfComponents="3" format="ascii">0.01 0.02 0.03</DataArray>
</PointData></Piece></PolyData></VTKFile>'''


class SnapshotAuditTests(unittest.TestCase):
    def scan(self, xml):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"LiquidCore_ite_0000000000.vtp"
            path.write_text(xml)
            return scan_vtp(path)

    def test_small_mass_survives_and_density_volume_reconciles(self):
        result = self.scan(VTP)
        self.assertEqual(result["points"], 1)
        self.assertEqual(result["mass_min_kg"], 9.350000000000002e-13)
        self.assertLess(result["density_volume_mass_relative_error"], 1e-14)

    def test_rounded_to_zero_particle_mass_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-positive particle Mass"):
            self.scan(VTP.replace("9.350000000000002e-13", "0.000000000"))

    def test_nonfinite_velocity_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "non-finite Velocity"):
            self.scan(VTP.replace("0.01 0.02 0.03", "0.01 NaN 0.03"))

    def test_declared_float32_is_rejected_in_precision_mode(self):
        with self.assertRaisesRegex(ValueError, "not Float64"):
            self.scan(VTP.replace('Name="Mass" type="Float64"', 'Name="Mass" type="Float32"'))

    def test_wrong_array_length_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "wrong tuple count"):
            self.scan(VTP.replace("0.01 0.02 0.03", "0.01 0.02"))

    def test_committed_nonfinite_csv_row_is_not_hidden(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/"progress.csv"
            p.write_text("time_s,value\n0,1\n1,NaN\n")
            with self.assertRaisesRegex(ValueError,"non-finite committed row"):
                rows(p)
            p.write_text("time_s,value\n0,1\n1,Na")
            self.assertEqual(rows(p),[{"time_s":0.,"value":1.}])


if __name__=="__main__":
    unittest.main()
