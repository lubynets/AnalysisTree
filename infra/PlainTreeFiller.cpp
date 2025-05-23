/* Copyright (C) 2019-2021 GSI, Universität Tübingen
   SPDX-License-Identifier: GPL-3.0-only
   Authors: Viktor Klochkov, Ilya Selyuzhenkov */
#include <algorithm>

#include "TTree.h"

#include "PlainTreeFiller.hpp"
#include "TaskManager.hpp"

namespace AnalysisTree {

void PlainTreeFiller::AddBranch(const std::string& branch_name) {
  branch_name_ = branch_name;
  in_branches_.emplace(branch_name);
}

void PlainTreeFiller::SetFieldsToIgnore(const std::vector<std::string>& fields_to_ignore) {
  if (branch_name_ == "") {
    throw std::runtime_error("PlainTreeFiller::SetFieldsToIgnore() must be called after PlainTreeFiller::AddBranch()\n");
  }
  for (auto& fti : fields_to_ignore) {
    fields_to_ignore_.emplace_back((branch_name_ + "." + fti).c_str());
  }
}

void PlainTreeFiller::SetFieldsToPreserve(const std::vector<std::string>& fields_to_preserve) {
  if (branch_name_ == "") {
    throw std::runtime_error("PlainTreeFiller::SetFieldsToPreserve() must be called after PlainTreeFiller::AddBranch()\n");
  }
  for (auto& fti : fields_to_preserve) {
    fields_to_preserve_.emplace_back((branch_name_ + "." + fti).c_str());
  }
}

void PlainTreeFiller::Init() {
  if (is_ignore_defual_fields_) {
    std::vector<std::string> defaultFieldsNames;
    auto mapF = config_->GetBranchConfig(branch_name_).GetMap<float>();
    auto mapI = config_->GetBranchConfig(branch_name_).GetMap<int>();
    auto mapB = config_->GetBranchConfig(branch_name_).GetMap<bool>();
    for (auto& m : {mapF, mapI, mapB}) {
      for (auto& me : m) {
        if (me.second.id_ < 0) defaultFieldsNames.emplace_back(me.first);
      }
    }
    SetFieldsToIgnore(std::move(defaultFieldsNames));
  }

  if (!fields_to_ignore_.empty() && !fields_to_preserve_.empty()) {
    throw std::runtime_error("PlainTreeFiller::Init() - only one of fields_to_ignore_ and fields_to_preserve_ can be set");
  }

  if (!branch_name_.empty()) {
    const auto& branch_config = config_->GetBranchConfig(branch_name_);
    for (const auto& field : branch_config.GetMap<float>()) {
      AnalysisTask::AddEntry(AnalysisEntry({Variable(branch_name_, field.first)}));
    }
    for (const auto& field : branch_config.GetMap<int>()) {
      AnalysisTask::AddEntry(AnalysisEntry({Variable(branch_name_, field.first)}));
    }
    for (const auto& field : branch_config.GetMap<bool>()) {
      AnalysisTask::AddEntry(AnalysisEntry({Variable(branch_name_, field.first)}));
    }
  }

  AnalysisTask::Init();

  if (entries_.size() != 1) {
    throw std::runtime_error("Only 1 output branch");
  }
  const auto& vars = entries_[0].GetVariables();
  vars_.resize(vars.size());

  file_ = TFile::Open(file_name_.c_str(), "recreate");
  plain_tree_ = new TTree(tree_name_.c_str(), "Plain Tree");
  plain_tree_->SetAutoSave(0);
  for (size_t i = 0; i < vars.size(); ++i) {
    std::string leaf_name = vars[i].GetName();
    if (!fields_to_ignore_.empty() && std::find(fields_to_ignore_.begin(), fields_to_ignore_.end(), leaf_name) != fields_to_ignore_.end()) continue;
    if (!fields_to_preserve_.empty() && std::find(fields_to_preserve_.begin(), fields_to_preserve_.end(), leaf_name) == fields_to_preserve_.end()) continue;
    if (!is_prepend_leaves_with_branchname_) leaf_name.erase(0, branch_name_.size() + 1);
    std::replace(leaf_name.begin(), leaf_name.end(), '.', '_');
    plain_tree_->Branch(leaf_name.c_str(), &(vars_.at(i)), Form("%s/F", leaf_name.c_str()));
  }

  for (auto& cm : cuts_map_) {
    if (cm.second != nullptr) {
      cm.second->Init(*(TaskManager::GetInstance()->GetChain()->GetConfiguration()));
    }
  }
}

void PlainTreeFiller::Exec() {
  AnalysisTask::Exec();
  const auto& values = entries_[0].GetValues();
  for (const auto& channel : values) {
    assert(channel.size() == vars_.size());
    for (size_t i = 0; i < channel.size(); ++i) {
      vars_.at(i) = channel.at(i);
    }
    plain_tree_->Fill();
  }
}

void PlainTreeFiller::Finish() {
  AnalysisTask::Finish();
  file_->cd();
  plain_tree_->Write();
  file_->Close();
  delete file_;
  //  delete plain_tree_;
}
}// namespace AnalysisTree
