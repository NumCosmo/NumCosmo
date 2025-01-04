#ifndef _NC_OBJECT_NAME_H_
#define _NC_OBJECT_NAME_H_

#include <glib.h>
#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_OBJECT_NAME (nc_object_name_get_type())

G_DECLARE_FINAL_TYPE(NcObjectName, nc_object_name, NC, OBJECT_NAME, GObject)

/* METHODS */
NcObjectName *nc_object_name_new(gdouble prop1_val);
NcObjectName *nc_object_name_ref(NcObjectName *nc_object_name);
void nc_object_name_free(NcObjectName *nc_object_name);
void nc_object_name_clear(NcObjectName **nc_object_name);

G_END_DECLS

#endif /* _NC_OBJECT_NAME_H_ */
