/*
 * mwcs.h
 *
 *  Created on: 1-feb-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCS_H
#define MWCS_H

typedef enum {
  MwcsSolverSCF = 0,
  MwcsSolverMCF,
  MwcsSolverCutFlow,
  MwcsSolverCutFlowMin,
  MwcsSolverCutNodeSeparator,
  MwcsSolverCutNodeSeparatorBk,
  MwcsSolverTreeDP,
  MwcsSizeSolverTreeDP,
  MwcsSizeSolverCutNodeSeparatorBk } MwcsSolverEnum;

#endif // MWCS_H
