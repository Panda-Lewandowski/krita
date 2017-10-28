/*
 * Copyright (C) 2016 Boudewijn Rempt <boud@valdyas.org>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this library; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef ColorModelPerLayerCHECK_H
#define ColorModelPerLayerCHECK_H

#include "KisExportCheckRegistry.h"
#include <KoID.h>
#include <klocalizedstring.h>
#include <kis_image.h>
#include <KoColorSpace.h>
#include <KoColorModelStandardIds.h>
#include <kis_layer.h>
#include <kis_node_visitor.h>
#include "kis_node.h"
#include "kis_paint_layer.h"
#include "kis_group_layer.h"
#include "kis_adjustment_layer.h"
#include "kis_external_layer_iface.h"
#include "kis_clone_layer.h"
#include "generator/kis_generator_layer.h"

class KisColorModelCheckVisitor : public KisNodeVisitor
{
public:

    using KisNodeVisitor::visit;

    KisColorModelCheckVisitor (KoID colorModelID, KoID colorDepthID)
        : m_count(0)
        , m_colorModelID(colorModelID)
        , m_colorDepthID(colorDepthID)
    {
    }

    quint32 count() {
        return m_count;
    }

    bool visit(KisNode* node) override {
        return check(node);
    }

    bool visit(KisPaintLayer *layer) override {
        return check(layer);
    }

    bool visit(KisGroupLayer *layer) override {
        return check(layer);
    }


    bool visit(KisAdjustmentLayer *layer) override {
        return check(layer);
    }

    bool visit(KisExternalLayer *layer) override {
        return check(layer);
    }

    bool visit(KisCloneLayer *layer) override {
        return check(layer);
    }

    bool visit(KisGeneratorLayer * layer) override {
        return check(layer);
    }

    bool visit(KisFilterMask *) override {return true;}

    bool visit(KisTransformMask *) override {return true;}

    bool visit(KisTransparencyMask *) override {return true;}

    bool visit(KisSelectionMask *) override {return true;}

    bool visit(KisColorizeMask *) override {return true;}


private:
    bool check(KisNode * node)
    {
        KisLayer *layer = dynamic_cast<KisLayer*>(node);
        if (layer) {
            const KoColorSpace * cs = layer->colorSpace();
            if (cs->colorModelId() == m_colorModelID && cs->colorDepthId() == m_colorDepthID) {
                m_count++;
            }
        }
        visitAll(node);
        return true;
    }

    quint32 m_count;
    const KoID m_colorModelID;
    const KoID m_colorDepthID;

};

class ColorModelPerLayerCheck : public KisExportCheckBase
{
public:

    ColorModelPerLayerCheck(const KoID &colorModelID, const KoID &colorDepthID, const QString &id, Level level, const QString &customWarning = QString())
        : KisExportCheckBase(id, level, customWarning, true)
        , m_ColorModelID(colorModelID)
        , m_colorDepthID(colorDepthID)
    {
        Q_ASSERT(!colorModelID.name().isEmpty());
        Q_ASSERT(!colorDepthID.name().isEmpty());

        if (customWarning.isEmpty()) {
            m_warning = i18nc("image conversion warning", "Your image contains layers with the color model <b>%1</b> and channel depth <b>%2</b> which cannot be saved to this format. The layers will be converted."
            , m_ColorModelID.name()
            ,m_colorDepthID.name());
        }
    }

    bool checkNeeded(KisImageSP image) const override
    {
        KisColorModelCheckVisitor v(m_ColorModelID, m_colorDepthID);
        image->rootLayer()->accept(v);
        return (v.count() > 0);
    }

    Level check(KisImageSP /*image*/) const override
    {
        return m_level;
    }

    const KoID m_ColorModelID;
    const KoID m_colorDepthID;
};

class ColorModelPerLayerCheckFactory : public KisExportCheckFactory
{
public:

    ColorModelPerLayerCheckFactory(const KoID &ColorModelID, const KoID &colorDepthId)
        : m_colorModelID(ColorModelID)
        , m_colorDepthID(colorDepthId)
    {
    }

    ~ColorModelPerLayerCheckFactory() override {}

    KisExportCheckBase *create(KisExportCheckBase::Level level, const QString &customWarning) override
    {
        return new ColorModelPerLayerCheck(m_colorModelID, m_colorDepthID, id(), level, customWarning);
    }

    QString id() const override {
        return "ColorModelPerLayerCheck/" + m_colorModelID.id() + "/" + m_colorDepthID.id();
    }

    const KoID m_colorModelID;
    const KoID m_colorDepthID;
};

#endif // ColorModelPerLayerCHECK_H
