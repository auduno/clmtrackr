import React, { PropTypes } from 'react';
import { connect } from 'react-redux';
import classNames from 'classnames';
import _ from 'lodash';

import Drawer from 'material-ui/Drawer';
import MenuItem from 'material-ui/MenuItem';

import {
  HEADER_TITLES,
  EXAMPLES,
  setExample
} from 'clmtrackr/examples/reducers/examples';

import './MenuDrawer.styl';


class MenuDrawer extends React.Component {
  render () {
    const sections = {};
    _.forEach(EXAMPLES, (example, i) => {
      const exampleType = example.type || 'example';
      if (!sections[exampleType]) {
        sections[exampleType] = [];
      }

      sections[exampleType].push(
        <MenuItem
          className={classNames({
            selected: example.id === this.props.selectedExample
          })}
          onClick={() => this.props.setExample(example)}
          key={example.id}
        >
          {example.name}
        </MenuItem>
      );
    });

    const menuItems = [];
    _.forEach(sections, (values, section) => {
      menuItems.push(
        <div className='section-header' key={section}>
          {HEADER_TITLES[section] || section}
        </div>
      );
      _.forEach(values, sectionEntry => menuItems.push(sectionEntry));
    })

    return (
      <div className='menu-drawer-cmpt'>
        <Drawer open={true}>
          {menuItems}
        </Drawer>
      </div>
    );
  }
}

MenuDrawer.propTypes = {
  selectedExample: PropTypes.string,
  setExample: PropTypes.func.isRequired
};

const mapStateToProps = state => {
  return {
    selectedExample: state.examples.activeExample
  };
};

const mapDispatchToProps = dispatch => {
  return {
    setExample: example => dispatch(setExample(example))
  };
};

export default connect(mapStateToProps, mapDispatchToProps)(MenuDrawer);
