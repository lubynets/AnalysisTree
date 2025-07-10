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
  if (branch_name_.empty()) {
    throw std::runtime_error("PlainTreeFiller::SetFieldsToIgnore() must be called after PlainTreeFiller::AddBranch()\n");
  }
  for (const auto& fti : fields_to_ignore) {
    fields_to_ignore_.emplace_back((branch_name_ + "." + fti).c_str());
  }
}

void PlainTreeFiller::SetFieldsToPreserve(const std::vector<std::string>& fields_to_preserve) {
  if (branch_name_.empty()) {
    throw std::runtime_error("PlainTreeFiller::SetFieldsToPreserve() must be called after PlainTreeFiller::AddBranch()\n");
  }
  for (const auto& fti : fields_to_preserve) {
    fields_to_preserve_.emplace_back((branch_name_ + "." + fti).c_str());
  }
}

void PlainTreeFiller::SetFieldsToRename(const std::vector<std::pair<std::string, std::string>>& fields_to_rename) {
  if (branch_name_.empty()) {
    throw std::runtime_error("PlainTreeFiller::SetFieldsToRename() must be called after PlainTreeFiller::AddBranch()\n");
  }
  for (const auto& ftr : fields_to_rename) {
    fields_to_rename_.emplace((branch_name_ + "." + ftr.first).c_str(), (branch_name_ + "." + ftr.second).c_str());
  }
}

void PlainTreeFiller::CheckIgnorePreserveRenameFields(const std::vector<std::string>& leafNames) const {
  for (const auto& fti : fields_to_ignore_) {
    if (std::find(leafNames.begin(), leafNames.end(), fti) == leafNames.end()) {
      std::cout << "WARNING PlainTreeFiller::CheckIgnorePreserveRenameFields(): field " << fti << " is set to be ignored, but it is absent among input fields\n";
    }
  }

  for (const auto& ftp : fields_to_preserve_) {
    if (std::find(leafNames.begin(), leafNames.end(), ftp) == leafNames.end()) {
      std::cout << "WARNING PlainTreeFiller::CheckIgnorePreserveRenameFields(): field " << ftp << " is set to be preserved, but it is absent among input fields\n";
    }
  }

  for (const auto& ftr : fields_to_rename_) {
    if (std::find(leafNames.begin(), leafNames.end(), ftr.first) == leafNames.end()) {
      std::cout << "WARNING PlainTreeFiller::CheckIgnorePreserveRenameFields(): field " << ftr.first << " is set to be renamed, but it is absent among input fields\n";
    }
  }
}

void PlainTreeFiller::Init() {
  if (is_ignore_defual_fields_) {
    std::vector<std::string> defaultFieldsNames;
    auto mapF = config_->GetBranchConfig(branch_name_).GetMap<float>();
    auto mapI = config_->GetBranchConfig(branch_name_).GetMap<int>();
    auto mapB = config_->GetBranchConfig(branch_name_).GetMap<bool>();
    for (const auto& m : {mapF, mapI, mapB}) {
      for (const auto& me : m) {
        if (me.second.id_ < 0) defaultFieldsNames.emplace_back(me.first);
      }
    }
    SetFieldsToIgnore(defaultFieldsNames);
  }

  if (!fields_to_ignore_.empty() && !fields_to_preserve_.empty()) {
    throw std::runtime_error("PlainTreeFiller::Init() - only one of fields_to_ignore_ and fields_to_preserve_ can be set");
  }

  if (!branch_name_.empty()) {
    const auto& branch_config = config_->GetBranchConfig(branch_name_);
    for (const auto& field : branch_config.GetMap<float>()) {
      AnalysisTask::AddEntry(AnalysisEntry({Variable(branch_name_, field.first)}));
      vars_.emplace_back();
      vars_.back().type_ = Types::kFloat;
    }
    for (const auto& field : branch_config.GetMap<int>()) {
      AnalysisTask::AddEntry(AnalysisEntry({Variable(branch_name_, field.first)}));
      vars_.emplace_back();
      vars_.back().type_ = Types::kInteger;
    }
    for (const auto& field : branch_config.GetMap<bool>()) {
      AnalysisTask::AddEntry(AnalysisEntry({Variable(branch_name_, field.first)}));
      vars_.emplace_back();
      vars_.back().type_ = Types::kBool;
    }
  }

  AnalysisTask::Init();

  if (entries_.size() != 1) {
    throw std::runtime_error("Only 1 output branch");
  }
  const auto& vars = entries_[0].GetVariables();

  if (vars_.size() != vars.size()) throw std::runtime_error("PlainTreeFiller::Init(): vars_.size() != vars.size()");

  std::vector<std::string> leafNames;
  for (int iVar = 0, nVars = vars.size(); iVar < nVars; ++iVar) {
    leafNames.emplace_back(vars[iVar].GetName());
  }
  CheckIgnorePreserveRenameFields(leafNames);

  file_ = TFile::Open(file_name_.c_str(), "recreate");
  plain_tree_ = new TTree(tree_name_.c_str(), "Plain Tree");
  plain_tree_->SetAutoSave(0);
  for (int iLeaf = 0, nLeafs = leafNames.size(); iLeaf < nLeafs; ++iLeaf) {
    std::string leaf_name = leafNames.at(iLeaf);
    if (!fields_to_ignore_.empty() && std::find(fields_to_ignore_.begin(), fields_to_ignore_.end(), leaf_name) != fields_to_ignore_.end()) continue;
    if (!fields_to_preserve_.empty() && std::find(fields_to_preserve_.begin(), fields_to_preserve_.end(), leaf_name) == fields_to_preserve_.end()) continue;
    if (!fields_to_rename_.empty() && fields_to_rename_.find(leaf_name) != fields_to_rename_.end()) {
      leaf_name = fields_to_rename_.at(leaf_name);
    }
    if (!is_prepend_leaves_with_branchname_) leaf_name.erase(0, branch_name_.size() + 1);
    std::replace(leaf_name.begin(), leaf_name.end(), '.', '_');
    if (vars_.at(iLeaf).type_ == Types::kFloat) plain_tree_->Branch(leaf_name.c_str(), &vars_.at(iLeaf).float_, Form("%s/F", leaf_name.c_str()));
    else if (vars_.at(iLeaf).type_ == Types::kInteger)
      plain_tree_->Branch(leaf_name.c_str(), &vars_.at(iLeaf).int_, Form("%s/I", leaf_name.c_str()));
    else if (vars_.at(iLeaf).type_ == Types::kBool)
      plain_tree_->Branch(leaf_name.c_str(), &vars_.at(iLeaf).bool_, Form("%s/O", leaf_name.c_str()));
  }

  for (const auto& cm : cuts_map_) {
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
      if (vars_.at(i).type_ == Types::kFloat) vars_.at(i).float_ = static_cast<float>(channel.at(i));
      else if (vars_.at(i).type_ == Types::kInteger)
        vars_.at(i).int_ = static_cast<int>(std::round(channel.at(i)));
      else if (vars_.at(i).type_ == Types::kBool)
        vars_.at(i).bool_ = static_cast<bool>(std::round(channel.at(i)));
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
}
}// namespace AnalysisTree
