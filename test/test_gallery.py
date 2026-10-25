"""Scientific and publication-boundary tests for the static gallery builder."""
import csv
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("build_gallery",ROOT/"scripts"/"build_gallery.py")
gallery=importlib.util.module_from_spec(spec)
spec.loader.exec_module(gallery)


def fixture(path, name="test-run", amplitude=1):
    path.mkdir()
    config={
        "experiment":{"name":name,"purpose":"maintenance","seed":42},
        "grid":{"nx":8,"length_m":8.0},
        "planet":{"radius_m":60.0,"rotation_rate_s":.000165},
        "physics":{"deformation_radius_m":2.0,"deep_jet_speed_m_s":.5,"deep_jet_radius_m":2.0,"deep_jet_width_m":1.0},
        "initial":{"jet_speed_m_s":1.0,"jet_radius_m":2.0,"jet_width_m":1.0,"cyclone_speed_m_s":.5,"cyclone_radius_m":1.0,"perturbation":"mode","displacement_m":.1,"mode":6,"m_max":6},
        "boundary":{"inner_radius_m":2.5,"outer_radius_m":3.0,"sponge_rate_s":0.0},
        "numerics":{"dt_s":.1,"duration_s":2.0,"hyper_order":4,"cutoff_damping_time_s":1e5,"cfl_limit":.6},
        "convection":{"enabled":False},"output":{"m_max":6}}
    (path/"config.toml").write_text(gallery.toml_text(config)+ '\n[private]\npath = "/Users/private/research"\n')
    (path/"metadata.toml").write_text('source_sha256 = "'+"a"*64+'"\nmanifest_sha256 = "'+"b"*64+'"\nresume_from="/Users/private/checkpoint.jld2"\nhostname="secret-node.example.org"\n')
    (path/"COMPLETED").write_text("2.0\n")
    x=np.arange(8,dtype=float)-4
    with h5py.File(str(path/"domain.jld2"),"w") as f:
        f["x"]=x;f["y"]=x
    k=2*math.pi/8
    psi=amplitude*np.sin(k*x[:,None])*np.cos(k*x[None,:])
    zeta=-2*k*k*psi
    fields={"psi":psi,"s":zeta-psi/4,"u":amplitude*k*np.sin(k*x[:,None])*np.sin(k*x[None,:]),"v":amplitude*k*np.cos(k*x[:,None])*np.cos(k*x[None,:]),"totalpv":zeta-psi/4}
    for i,time in enumerate([0.0,2.0]):
        with h5py.File(str(path/("snapshot_%09d.jld2" % i)),"w") as f:
            f["time_s"]=time
            for key,value in fields.items():f[key]=value.T
    keys=gallery.CORE_COLUMNS+["m%d_%s_m"%(m,part) for m in range(1,7) for part in ("amplitude","real","imag")]+["private_path"]
    with (path/"diagnostics.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(keys)
        for time in [0,1,2]:
            row=[time,time*10,1,1,0,2,1,1]
            for m in range(1,7):row.extend([.1*m,.05*m,0])
            writer.writerow(row+["/Users/private/diagnostic"])
    return config,fields,zeta


class GalleryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.root=Path(self.tmp.name)
        self.result=self.root/"result"
        self.config,self.fields,self.zeta=fixture(self.result)
    def tearDown(self):self.tmp.cleanup()

    def test_analytic_vorticity_and_barotropic_limit(self):
        np.testing.assert_allclose(gallery.vorticity(self.fields,2),self.zeta,atol=1e-14)
        np.testing.assert_array_equal(gallery.vorticity(self.fields,math.inf),self.fields["s"])
        with self.assertRaises(ValueError):gallery.vorticity(self.fields,0)

    def test_whole_deep_profile_is_retained_and_labelled(self):
        self.config["physics"].update(deep_profile="matched_upper", deep_match_fraction=.5, deep_jet_speed_m_s=0)
        public = gallery.safe_config(self.config)
        self.assertEqual(public["physics"]["deep_match_fraction"], .5)
        self.assertIn("jet + cyclone", gallery.deep_flow_label(public["physics"]))
        self.assertIn("0.5", gallery.deep_flow_label(public["physics"]))
        self.assertEqual(gallery.deep_flow_label({"deep_jet_speed_m_s":50}), "Annular jet 50 m/s")
        self.config["physics"]["deep_profile"] = "/private/path"
        with self.assertRaises(ValueError):gallery.safe_config(self.config)

    def test_hdf5_and_display_orientation(self):
        time,fields=gallery.read_snapshot(self.result/"snapshot_000000000.jld2",(8,8))
        np.testing.assert_array_equal(fields["psi"],self.fields["psi"])
        asymmetric=np.arange(8)[:,None]+10*np.arange(8)[None,:]
        image,factor=gallery.display_array(asymmetric,8)
        self.assertEqual(image[2,5],asymmetric[5,2])
        self.assertEqual(factor,1)
        reduced,factor=gallery.display_array(asymmetric,4)
        self.assertEqual(reduced[1,2],asymmetric[4:6,2:4].mean())

    def test_times_reject_nonfinite_duplicates_and_reversal(self):
        for times in ([0,np.nan],[0,np.inf],[0,0],[2,1],[-1,0],[]):
            with self.assertRaises(ValueError):gallery.increasing_times(times,"test")
        with h5py.File(str(self.result/"snapshot_000000001.jld2"),"a") as f:f["time_s"][...]=0
        with self.assertRaises(ValueError):gallery.load_run(self.result)

    def test_nonfinite_field_rejected(self):
        with h5py.File(str(self.result/"snapshot_000000001.jld2"),"a") as f:f["u"][0,0]=np.nan
        with self.assertRaises(ValueError):gallery.load_run(self.result)

    def test_completion_and_hash_validation(self):
        (self.result/"COMPLETED").write_text("1")
        with self.assertRaises(ValueError):gallery.load_run(self.result)
        (self.result/"COMPLETED").write_text("2")
        (self.result/"Manifest.toml").write_text("not matching recorded hash")
        with self.assertRaises(ValueError):gallery.load_run(self.result)

    def test_refuses_overwrite_before_loading(self):
        output=self.root/"existing";output.mkdir();(output/"keep").write_text("untouched")
        with self.assertRaises(FileExistsError):gallery.build_gallery([self.root/"absent"],output)
        self.assertEqual((output/"keep").read_text(),"untouched")

    def test_missing_modes_remain_gaps(self):
        path=self.result/"diagnostics.csv"
        with path.open() as stream:
            rows=list(csv.DictReader(stream))
        rows[1]["jet_coverage"]="0.5";rows[1]["m6_amplitude_m"]="nan"
        with path.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        data=gallery.read_diagnostics(path)
        self.assertTrue(np.isnan(data["m6_amplitude_m"][1]))
        self.assertNotIn("private_path",data)

    def test_broadband_and_absent_disturbance_labels(self):
        run=gallery.load_run(self.result)
        initial=run["public_config"]["initial"]
        initial.update(perturbation="broadband",displacement_m=10000,m_max=32)
        rows=dict(gallery.metadata_rows(run))
        self.assertIn("normalization d = 10 km",rows["Initial disturbance"])
        self.assertIn("per-mode parameter d/√M = 1.7678 km",rows["Initial disturbance"])
        self.assertEqual(rows["Initial azimuthal mode(s)"],"1–32")
        initial["perturbation"]="none"
        rows=dict(gallery.metadata_rows(run))
        self.assertEqual(rows["Initial disturbance"],"None")
        self.assertEqual(rows["Initial azimuthal mode(s)"],"None")

    def test_plot_title_when_mode_six_is_not_recorded(self):
        run=gallery.load_run(self.result)
        run["data"]={key:value for key,value in run["data"].items()
                     if not gallery.MODE_COLUMN.fullmatch(key) or int(gallery.MODE_COLUMN.fullmatch(key).group(1))<=2}
        with patch.object(gallery.plt,"close"):
            gallery.render_diagnostics(run,self.root/"two_modes.png")
            figure=gallery.plt.gcf()
            self.assertEqual(figure.axes[0].get_title(),"Measured mode amplitudes")
        gallery.plt.close(figure)

    def test_actual_build_escaping_shared_limits_and_privacy(self):
        second=self.root/"second";fixture(second,name="/Users/private/run",amplitude=2)
        output=self.root/"gallery"
        title='<script>alert("bad")</script> & science'
        manifest=gallery.build_gallery([self.result,second],output,title=title)
        self.assertEqual(len(manifest["runs"]),2)
        self.assertAlmostEqual(manifest["shared_color_limits"]["speed_m_s"][1],math.pi/2)
        index=(output/"index.html").read_text()
        self.assertNotIn('<script>alert',index)
        self.assertIn('&lt;script&gt;',index)
        self.assertTrue((output/"COMPLETED").exists())
        for file in output.rglob("*"):
            if file.is_file() and file.suffix in (".html",".json",".toml",".csv",".js",".css"):
                text=file.read_text()
                for private in ("/Users/private",str(self.root),"secret-node.example.org"):
                    self.assertNotIn(private,text)
        for run in manifest["runs"]:
            folder=output/"runs"/run["run_id"]
            summary=json.loads((folder/"summary.json").read_text())
            self.assertEqual([f["time_s"] for f in summary["frames"]],[0,2])
            self.assertEqual(summary["shared_color_limits"],manifest["shared_color_limits"])
            self.assertEqual(gallery.sha256(folder/"config.toml"),summary["public_config_sha256"])
            self.assertEqual(len(list(folder.glob("frame_*.png"))),2)
            self.assertFalse(list(folder.glob("*.jld2")))


if __name__=="__main__":unittest.main()
