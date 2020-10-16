//-----------------------------------------------------------------------------
//   qualrout.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             04/02/15   (Build 5.1.008)
//             04/30/15   (Build 5.1.009)
//             08/05/15   (Build 5.1.010)
//             04/01/20   (Build 5.1.015)
//   Author:   L. Rossman
//
//   Water quality routing functions.
//
//   Build 5.1.008:
//   - Pollutant mass lost to seepage flow added to mass balance totals.
//   - Pollutant concen. increased when evaporation occurs.
//
//   Build 5.1.009:
//   - Criterion for dry link/storage node changed to avoid concen. blowup.
//
//   Build 5.1.010:
//   - Entire module re-written to be more compact and easier to follow.
//   - Neglible depth limit replaced with a negligible volume limit.
//
//   Build 5.1.015:
//   - Fixed mass balance issue for empty storage nodes that flood.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const double ZeroVolume = 0.0353147; // 1 liter in ft3

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  qualrout_init            (called by routing_open)
//  qualrout_execute         (called by routing_execute)

//-----------------------------------------------------------------------------
//  Function declarations
//-----------------------------------------------------------------------------
static void  findLinkMassFlow(int i, double tStep);
static void  findNodeQual(int j);
static void  findLinkQual(int i, double tStep);
static void  findSFLinkQual(int i, double qSeep, double fEvap, double tStep);
static void  findStorageQual(int j, double tStep);
static void  updateHRT(int j, double v, double q, double tStep);
static double getReactedQual(int p, double c, double v1, double tStep);
static double getMixedQual(double c, double v1, double wIn, double qIn,
              double tStep);
/* START modification by Peter Schlagbauer | TUGraz */
//static double getReactedTemp(double oldTemp, int i, double tStep);
//static double getReactedTempStNode(double oldTemp, int j, int p, double tStep);
//static double getWettedArea(TTable* table, double d);
/* END modification by Peter Schlagbauer | TUGraz */
//=============================================================================

void    qualrout_init()
//
//  Input:   none
//  Output:  none
//  Purpose: initializes water quality concentrations in all nodes and links.
//
{
    int     i, p, isWet;
    double  c;

    for (i = 0; i < Nobjects[NODE]; i++)
    {
        isWet = ( Node[i].newDepth > FUDGE );
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
        /* START modification by Peter Schlagbauer | TUGraz */
        //if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1)
		//c = NAN; // set temperature to NaN, because 0 is a valid temperatur value
	//else
		c = 0.0;
	/* END modification by Peter Schlagbauer | TUGraz */
            if ( isWet ) c = Pollut[p].initConcen;
            Node[i].oldQual[p] = c;
            Node[i].newQual[p] = c;
        }
	//if (Node[i].type == STORAGE)
	//	{
	//		int k = Node[i].subIndex;
			/* START modification by Peter Schlagbauer | TUGraz */
			// calculate the penetration depth for each storage node
		//	Storage[k].penDepth = sqrt(Storage[k].kSoil / (0.0000727220 * (Storage[k].densitySoil * Storage[k].specHcSoil)));
			/* END modification by Peter Schlagbauer | TUGraz */
	//	}
    }

    for (i = 0; i < Nobjects[LINK]; i++)
    {
        isWet = ( Link[i].newDepth > FUDGE );
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
		/* START modification by Peter Schlagbauer | TUGraz */
		//if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1)
		//	c = NAN; // set temperature to NaN, because 0 is a valid temperatur value
		//else
			c = 0.0;
		/* END modification by Peter Schlagbauer | TUGraz */
            if ( isWet ) c = Pollut[p].initConcen;
            Link[i].oldQual[p] = c;
            Link[i].newQual[p] = c;
        }

	//int k = Link[i].subIndex;
	/* START modification by Peter Schlagbauer | TUGraz */
	// calculate the penetration depth for each conduit
	//Conduit[k].penDepth = sqrt(Conduit[k].kSoil / (0.0000727220 * Conduit[k].densitySoil * Conduit[k].specHcSoil));
    //fprintf(stdout,"%g\n", Conduit[k].penDepth);
    /* END modification by Peter Schlagbauer | TUGraz */
    }
}

//=============================================================================

void qualrout_execute(double tStep)
//
//  Input:   tStep = routing time step (sec)
//  Output:  none
//  Purpose: routes water quality constituents through the drainage
//           network over the current time step.
//
{
    int    i, j;
    double qIn, vAvg;
    double unitScale;

    // get the current month of simulation
    //DateTime currentDate = getDateTime(NewRoutingTime);
    //int month = datetime_monthOfYear(currentDate);


    // --- find mass flow each link contributes to its downstream node
    for ( i = 0; i < Nobjects[LINK]; i++ ) findLinkMassFlow(i, tStep);

    // --- find new water quality concentration at each node  
    for (j = 0; j < Nobjects[NODE]; j++)
    {
        // --- get node inflow and average volume
        qIn = Node[j].inflow;
        vAvg = (Node[j].oldVolume + Node[j].newVolume) / 2.0;
        
        // --- save inflow concentrations if treatment applied
        if ( Node[j].treatment )
        {
            if ( qIn < ZERO ) qIn = 0.0;
            treatmnt_setInflow(qIn, Node[j].newQual);
        }
       
        // --- find new quality at the node 
        if ( Node[j].type == STORAGE || Node[j].oldVolume > FUDGE )
        {
            //findStorageQual(j, tStep, month);
            findStorageQual(j, tStep);
        }
        else findNodeQual(j);

        // --- apply treatment to new quality values
        if ( Node[j].treatment ) treatmnt_treat(j, qIn, vAvg, tStep);
    }

    // --- find new water quality in each link
    //unitScale = UCF(LENGTH) / (2.0 * UCF(FLOW));
    //for ( i = 0; i < Nobjects[LINK]; i++ ) findLinkQual(i, tStep, unitScale, month);
    for (i = 0; i < Nobjects[LINK]; i++) findLinkQual(i, tStep);
}

//=============================================================================

double getMixedQual(double c, double v1, double wIn, double qIn, double tStep)
//
//  Input:   c = concentration in reactor at start of time step (mass/ft3)
//           v1 = volume in reactor at start of time step (ft3)
//           wIn = mass inflow rate (mass/sec)
//           qIn = flow inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  returns pollutant concentration at end of time step (mass/ft3)
//  Purpose: finds pollutant concentration within a completely mixed reactor.
//
{
    double vIn, cIn, cMax;

    // --- if no inflow then reactor concentration is unchanged
    if ( qIn <= ZERO ) return c;

    // --- compute concentration of any inflow
    vIn = qIn * tStep;
    cIn = wIn * tStep / vIn;

    // --- mixture concen. can't exceed either original or inflow concen.
    cMax = MAX(c, cIn);

    // --- mix inflow with current reactor contents
    c = (c*v1 + wIn*tStep) / (v1 + vIn);
    c = MIN(c, cMax);
    c = MAX(c, 0.0);
    return c;
}


//=============================================================================

void findLinkMassFlow(int i, double tStep)
//
//  Input:   i = link index
//           tStep = time step (sec)
//  Output:  none
//  Purpose: adds constituent mass flow out of link to the total
//           accumulation at the link's downstream node.
//
//  Note:    Node[].newQual[], the accumulator variable, already contains
//           contributions from runoff and other external inflows from
//           calculations made in routing_execute().
{
    int    j, p;
    double qLink, w;

    // --- find inflow to downstream node
    qLink = Link[i].newFlow;

    // --- identify index of downstream node
    j = Link[i].node2;
    if ( qLink < 0.0 ) j = Link[i].node1;
    qLink = fabs(qLink);

    // --- examine each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
	/* START modification by Peter Schlagbauer | TUGraz */
        // --- temporarily accumulate inflow load in Node[j].newQual
       // if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1)
       // {
       //     if (!isnan(Link[i].oldQual[p])) // do not consider NaN values
       //     {
                w = qLink * Link[i].oldQual[p];
                Node[j].newQual[p] += w;
        //    }
        //    else
        //        w = 0;
       // }
       // else
       // {
        //    w = qLink * Link[i].oldQual[p];
       //     Node[j].newQual[p] += w;
       // }
	/* END modification by Peter Schlagbauer | TUGraz */

        // --- update total load transported by link
        Link[i].totalLoad[p] += w * tStep;
    }
}

//=============================================================================

void findNodeQual(int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: finds new quality in a node with no storage volume.
//
{
    int    p;
    double qNode;

    // --- if there is flow into node then concen. = mass inflow/node flow
    qNode = Node[j].inflow;
    if ( qNode > ZERO )
    {
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
            Node[j].newQual[p] /= qNode;
        }
    }
    
    // --- otherwise concen. is 0
    else for (p = 0; p < Nobjects[POLLUT]; p++){
		/* START modification by Peter Schlagbauer | TUGraz */
		// set temperature/newQual to NaN, because 0 is a valid temperatur value
	//	if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1)
	//		Node[j].newQual[p] = NAN;
    //    else
            Node[j].newQual[p] = 0.0;
		/* END modification by Peter Schlagbauer | TUGraz */
	}
}

//=============================================================================

void findLinkQual(int i, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a link at end of the current time step.
//
{
    int    j,                // upstream node index
           k,                // conduit index
           p;                // pollutant index
    double wIn,              // pollutant mass inflow rate (mass/sec)
           qIn,              // inflow rate (cfs)
           qSeep,            // rate of seepage loss (cfs)
           v1,               // link volume at start of time step (ft3)
           v2,               // link volume at end of time step (ft3)
           c1,               // current concentration within link (mass/ft3)
           c2,               // new concentration within link (mass/ft3)
           vEvap,            // volume lost to evaporation (ft3)
           vLosses,          // evap. + seepage volume loss (ft3)
           fEvap,            // evaporation concentration factor
           barrels;          // number of barrels in conduit

    // --- identify index of upstream node
    j = Link[i].node1;
    if ( Link[i].newFlow < 0.0 ) j = Link[i].node2;

    // --- link quality is that of upstream node when
    //     link is not a conduit or is a dummy link
    if ( Link[i].type != CONDUIT || Link[i].xsect.type == DUMMY )
    {
        for (p = 0; p < Nobjects[POLLUT]; p++)
        {
            Link[i].newQual[p] = Node[j].newQual[p];
        }
        return;
    }

    // --- get flow rates and evaporation loss
    k = Link[i].subIndex;
    barrels = Conduit[k].barrels;
    qIn  = fabs(Conduit[k].q1) * barrels;
    qSeep = Conduit[k].seepLossRate * barrels;
    vEvap = Conduit[k].evapLossRate * barrels * tStep;

    // --- get starting and ending volumes
    v1 = Link[i].oldVolume;
    v2 = Link[i].newVolume;
    vLosses = qSeep*tStep + vEvap;

    // --- compute factor by which concentrations are increased due to
    //     evaporation loss 
    fEvap = 1.0;
    if ( vEvap > 0.0 && v1 > ZeroVolume ) fEvap += vEvap / v1;

    // --- Steady Flow routing requires special treatment
    if ( RouteModel == SF )
    {
        findSFLinkQual(i, qSeep, fEvap, tStep);
        return;
    }

    // --- adjust inflow to compensate for volume change under Dynamic
    //     Wave routing (which produces just a single (out)flow rate
    //     for a conduit)
    if ( RouteModel == DW )
    {
        qIn = qIn + (v2 + vLosses - v1) / tStep; 
        qIn = MAX(qIn, 0.0);
    }

    // --- examine each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
        // --- start with concen. at start of time step
        c1 = Link[i].oldQual[p];

        // --- update mass balance accounting for seepage loss
        massbal_addSeepageLoss(p, qSeep*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;
	/* START modification by Peter Schlagbauer | TUGraz */
		// it has been observed that at low flow rates the model may become unstable, therefore 0.5 L/s is a boundary
		//if (Link[i].newFlow * UCF(FLOW) / 1000 > 0.0005) 
	//if (Link[i].newFlow > unitScale * Link[i].xsect.yFull)
		//if (Link[i].newFlow * UCF(FLOW)  > 0.0005) 
	//	{
	//		if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1)
	//		{
				// --- adjust temperature by heat exchange processes
	//			if (Node[j].newQual[p] > 0.0) 					// one ore more inflows into the node
	//				c2 = getReactedTemp(Node[j].newQual[p], i, tStep, month);
	//			else // no inflow into the node, but still water inside the conduit
	//				c2 = getReactedTemp(c1, i, tStep, month);
	//		}
	//		else
	//		{
				// --- reduce concen. by 1st-order reaction
				c2 = getReactedQual(p, c1, v1, tStep);
	//		}
	//	}
	//	else
	//		c2 = c1;

		// --- mix resulting contents with inflow from upstream node
		//if (!isnan(Node[j].newQual[p])) { // do not consider NaN values
			wIn = Node[j].newQual[p] * qIn;
			c2 = getMixedQual(c2, v1, wIn, qIn, tStep);
		//}

		// --- set concen. to zero if remaining volume is negligible
		if (v2 < ZeroVolume)
		{
			massbal_addToFinalStorage(p, c2 * v2);

			// set temperature to NaN, because 0 is a valid temperatur value
			//if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1)
			//	c2 = NAN;
			//else
				c2 = 0.0;
		}
		/* END modification by Peter Schlagbauer | TUGraz */

     	        // --- assign new concen. to link
      	        Link[i].newQual[p] = c2;
    }
}

//=============================================================================

void  findSFLinkQual(int i, double qSeep, double fEvap, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a link at end of the current time step for
//           Steady Flow routing.
//
{
    int j = Link[i].node1;
    int p;
    double c1, c2;
    double lossRate;

    // --- examine each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
        // --- conduit's quality equals upstream node quality
        c1 = Node[j].newQual[p];

        // --- update mass balance accounting for seepage loss
        massbal_addSeepageLoss(p, qSeep*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- apply first-order decay over travel time
        c2 = c1;
        if ( Pollut[p].kDecay > 0.0 )
        {
            c2 = c1 * exp(-Pollut[p].kDecay * tStep);
            c2 = MAX(0.0, c2);
            lossRate = (c1 - c2) * Link[i].newFlow;
            massbal_addReactedMass(p, lossRate);
        }
        Link[i].newQual[p] = c2;
    }
}

//=============================================================================

void  findStorageQual(int j, double tStep)
//
//  Input:   j = node index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a node with storage volume.
//  
{
    int    p,                // pollutant index
           k;                // storage unit index
    double qIn,              // inflow rate (cfs)
           wIn,              // pollutant mass inflow rate (mass)
           v1,               // volume at start of time step (ft3)
           c1,               // initial pollutant concentration (mass/ft3)
           c2,               // final pollutant concentration (mass/ft3)
           qExfil = 0.0,     // exfiltration rate from storage unit (cfs)
           vEvap = 0.0,      // evaporation loss from storage unit (ft3)
           fEvap = 1.0;      // evaporation concentration factor

    // --- get inflow rate & initial volume
    qIn = Node[j].inflow;
    v1 = Node[j].oldVolume;

    // -- for storage nodes
    if ( Node[j].type == STORAGE )
    {    
        // --- update hydraulic residence time
        //     (HRT can be used in treatment functions)
        updateHRT(j, Node[j].oldVolume, qIn, tStep);

        // --- get exfiltration rate and evaporation loss
        k = Node[j].subIndex;
        qExfil = Storage[k].exfilLoss / tStep;
        vEvap = Storage[k].evapLoss;

        // --- compute factor by which concentrations are increased due to
        //     evaporation loss (avoiding huge factors as storage unit
        //     dries out completely)
        if ( vEvap > 0.0 && v1 > ZeroVolume ) fEvap += vEvap / v1;
    }

    // --- for each pollutant
    for (p = 0; p < Nobjects[POLLUT]; p++)
    {
        // --- start with concen. at start of time step 
        c1 = Node[j].oldQual[p];

        // --- update mass balance accounting for exfiltration loss
        massbal_addSeepageLoss(p, qExfil*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;
		/* START modification by Peter Schlagbauer | TUGraz */
		//if (strcmp(QualUnitsWords[Pollut[p].units], "CELSIUS") == 0 && TempModel.active == 1 && c1 != 0.0 && !isnan(c1))
		//	c1 = getReactedTempStNode(c1, j, p, tStep, month);
		//else
		//{
			// --- apply first order reaction only if no separate treatment function
			if (Node[j].treatment == NULL ||
				Node[j].treatment[p].equation == NULL)
			{
				c1 = getReactedQual(p, c1, v1, tStep);
			}
		//}
		/* END modification by Peter Schlagbauer | TUGraz */

        // --- mix resulting contents with inflow from all sources
        //     (temporarily accumulated in Node[j].newQual)
        wIn = Node[j].newQual[p];
        c2 = getMixedQual(c1, v1, wIn, qIn, tStep);

// --- set concen. to zero if remaining volume & inflow is negligible          //(5.1.015)
        if (Node[j].newVolume <= ZeroVolume && qIn <= FLOW_TOL)                //(5.1.015)
        {
            massbal_addToFinalStorage(p, c2 * Node[j].newVolume);
            c2 = 0.0;
        }

        // --- assign new concen. to node
        Node[j].newQual[p] = c2;
    }
}

//=============================================================================

void updateHRT(int j, double v, double q, double tStep)
//
//  Input:   j = node index
//           v = storage volume (ft3)
//           q = inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: updates hydraulic residence time (i.e., water age) at a 
//           storage node.
//
{
    int    k = Node[j].subIndex;
    double hrt = Storage[k].hrt;
    if ( v < ZERO ) hrt = 0.0;
    else hrt = (hrt + tStep) * v / (v + q*tStep);
    Storage[k].hrt = MAX(hrt, 0.0);
}

//=============================================================================

double getReactedQual(int p, double c, double v1, double tStep)
//
//  Input:   p = pollutant index
//           c = initial concentration (mass/ft3)
//           v1 = initial volume (ft3)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: applies a first order reaction to a pollutant over a given
//           time step.
//
{
    double c2, lossRate;
    double kDecay = Pollut[p].kDecay;

    if ( kDecay == 0.0 ) return c;
    c2 = c * (1.0 - kDecay * tStep);
    c2 = MAX(0.0, c2);
    lossRate = (c - c2) * v1 / tStep;
    massbal_addReactedMass(p, lossRate);
    return c2;
}
