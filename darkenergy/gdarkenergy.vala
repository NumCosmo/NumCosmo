/***************************************************************************
 *            gdarkenergy.vala
 *
 *  Fri December 21 15:11:00 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

using GLib;
using Gtk;

public class Nc.GDE : GLib.Object {
  private Gtk.Builder builder;
  private Gtk.Window main_window;

  [CCode (instance_pos = -1)]
  public void on_main_window_destroy (Gtk.Window source) {
    Gtk.main_quit ();
  }
  
  public GDE () {
      builder = new Builder ();
  }

  ~GDE () {
    builder = null;
    main_window = null;
  }

  public bool load_ui () {
    try {
      builder.add_from_file ("/home/sandro/Projects/numcosmo/darkenergy/gdarkenergy.glade");
      builder.connect_signals (this);
      main_window = builder.get_object ("main_window") as Gtk.Window;
      main_window.show_all ();
      return true;
    } catch (Error e) {
      stderr.printf ("Could not load UI: %s\n", e.message);
      return false;
    }
  }

  public void add_model_by_type (Gtk.ListStore store, Type pt) {
    foreach (unowned GLib.Type ch in pt.children ()) {
      if (ch.is_abstract ()) {
        add_model_by_type (store, ch);
      }
      else {
        Gtk.TreeIter iter;
        store.append (out iter);
        store.set (iter, 0, ch.name ());
      }
    }    
  }
  
  public void load_models () {
    Type mt = typeof (Nc.HICosmo);
    var model_list = builder.get_object ("model_list") as Gtk.ListStore;
    add_model_by_type (model_list, mt);
  }
  
  static int main (string[] args) {
    Gtk.init (ref args);
    Ncm.cfg_init ();
    var app = new Nc.GDE ();
    if (!app.load_ui ())
      return 1;

    app.load_models ();
    
    Gtk.main ();
    return 0;
  }
}
