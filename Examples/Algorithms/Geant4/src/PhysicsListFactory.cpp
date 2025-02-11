// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/PhysicsListFactory.hpp"

#include <G4VUserPhysicsList.hh>

namespace ActsExamples::Geant4 {

PhysicsListFactoryFunction::PhysicsListFactoryFunction(Function function)
    : m_function(std::move(function)) {}

std::unique_ptr<G4VUserPhysicsList> PhysicsListFactoryFunction::factorize()
    const {
  return m_function();
}

}  // namespace ActsExamples::Geant4
